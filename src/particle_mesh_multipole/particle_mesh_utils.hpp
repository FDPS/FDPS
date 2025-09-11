#pragma once
#include "../ps_defs.hpp"
//#include "multidimensional_array.hpp"
#include "fmm.hpp"
#include "cell.hpp"

namespace ParticleSimulator
{

    template <int p>
    class GreenFunction
    {
    public:
        typedef double real_t;
        // typedef double _Complex      cplx_t;
        typedef std::complex<real_t> cplx_t;
        // [TODO]
        //    GreenFunction_PBC in the original PMMM uses std::complex<real_t>
        //    as cplx_t. What is the difference between double _Complex and
        //    std::complex<double> ?
        //    --> When using double _Complex, a compilation error occurs
        //        the following code:
        //        const cplx_t phase(cos(theta), sin(theta));
        enum
        {
            LEN = lmbuf<p>::length,
            LEN2 = lmbuf<2 * p>::length,
        };

        int bc;
        int NX, NY, NZ;

        MultidimensionalArray<real_t, 4> gf_r;
        MultidimensionalArray<cplx_t, 4> gf_k;

    private:
        inline void initOpenImpl()
        {
            gf_r.initialize(2 * NZ, 2 * NY, 2 * NX, LEN2);
            gf_k.initialize(2 * NZ, 2 * NY, 1 + NX, LEN2);
        }

        inline void initPeriodicXYZImpl()
        {
            gf_r.initialize(NZ, NY, NX, LEN2);
            gf_k.initialize(NZ, NY, 1 + NX / 2, LEN2);
        }

        void setOpenImpl(const int icut, const double cell_length)
        {
            int i, j, k;
            for (k = 0; k < 2 * NZ; k++)
                for (j = 0; j < 2 * NY; j++)
                    for (i = 0; i < 2 * NX; i++)
                    {
                        if (k == NZ || j == NY || i == NX)
                        {
                            for (int lm = 0; lm < LEN2; lm++)
                            {
                                gf_r(k, j, i, lm) = 0.0;
                            }
                            continue;
                        }
                        const int kk = (k > NZ) ? k - 2 * NZ : k;
                        const int jj = (j > NY) ? j - 2 * NY : j;
                        const int ii = (i > NX) ? i - 2 * NX : i;
                        if (abs(kk) <= icut && abs(jj) <= icut && abs(ii) <= icut)
                        {
                            for (int lm = 0; lm < LEN2; lm++)
                            {
                                gf_r(k, j, i, lm) = 0.0;
                            }
                            continue;
                        }
                        const double dx = cell_length * double(ii);
                        const double dy = cell_length * double(jj);
                        const double dz = cell_length * double(kk);
                        Slm<2 * p, real_t> slm;
                        slm.eval_opt(-dx, dy, dz); // eval S_l^{-m}
                        for (int lm = 0; lm < LEN2; lm++)
                        {
                            gf_r(k, j, i, lm) = slm.buf[lm];
                        }
                    }
        }

        void setPeriodicXYZImpl(const int icut,
                                const double cell_length,
                                const double alpha,
                                const int NMAX,
                                const int MMAX)
        {
            {
                Slm<2 * p, real_t> dry;
                dry.init();
            }
#pragma omp parallel for
            for (int k = 0; k < NZ; k++)
                for (int j = 0; j < NY; j++)
                    for (int i = 0; i < NX; i++)
                    {
                        CutFunc<2 * p> cf;
                        // real-space sum
                        Slm<2 * p, real_t> rsum;
                        rsum.clear();

                        const int kk = (k >= NZ / 2) ? k - NZ : k;
                        const int jj = (j >= NY / 2) ? j - NY : j;
                        const int ii = (i >= NX / 2) ? i - NX : i;
                        int nx, ny, nz;
                        for (nz = -NMAX; nz <= NMAX; nz++)
                            for (ny = -NMAX; ny <= NMAX; ny++)
                                for (nx = -NMAX; nx <= NMAX; nx++)
                                {
                                    const int kkk = kk + nz * NZ;
                                    const int jjj = jj + ny * NY;
                                    const int iii = ii + nx * NX;
                                    if (0 == (iii | jjj | kkk))
                                        continue;
                                    const double dx = cell_length * double(iii);
                                    const double dy = cell_length * double(jjj);
                                    const double dz = cell_length * double(kkk);
                                    const double dr2 = dx * dx + dy * dy + dz * dz;
                                    cf.eval_rcut(dr2 * (alpha * alpha));
                                    Slm<2 * p, real_t> slm;
                                    slm.eval_opt(-dx, dy, dz); // eval S_l^{-m}
                                    // near cell correction
                                    const bool near = (abs(kkk) <= icut && abs(jjj) <= icut && abs(iii) <= icut);
                                    if (near)
                                    {
                                        for (int l = 0; l <= 2 * p; l++)
                                        {
                                            cf.rcut[l] -= 1.0;
                                        }
                                    }
                                    for (int l = 0; l <= 2 * p; l++)
                                    {
                                        for (int m = 0; m <= l; m++)
                                        {
                                            const cplx_t val = cf.rcut[l] * slm.val_at(l, m);
                                            rsum.accum_at(l, m, val);
                                        }
                                    }
                                }

                        // wave-space sum
                        Slm<2 * p, real_t> ksum;
                        ksum.clear();

                        int mx, my, mz;
                        for (mz = -MMAX; mz <= MMAX; mz++)
                            for (my = -MMAX; my <= MMAX; my++)
                                for (mx = -MMAX; mx <= MMAX; mx++)
                                {
                                    if (0 == (mx | my | mz))
                                        continue;
                                    const double dx = cell_length * double(i);
                                    const double dy = cell_length * double(j);
                                    const double dz = cell_length * double(k);

                                    const double theta = (+8.0 * atan(1.0)) * (dx * mx + dy * my + dz * mz); // 2pi * k \cdot r
                                    const cplx_t phase(cos(theta), sin(theta));

                                    Slm<2 * p, real_t> slm;
                                    slm.eval_opt(-double(mx), double(my), double(mz));

                                    const double m2 = mx * mx + my * my + mz * mz;
                                    cf.eval_kcut(m2, alpha);

                                    for (int l = 0; l <= 2 * p; l++)
                                    {
                                        for (int m = 0; m <= l; m++)
                                        {
                                            const cplx_t val = (cf.kcut[l] * phase) * slm.val_at(l, m);
                                            ksum.accum_at(l, m, val);
                                        }
                                    }
                                }
                                // store sum
#ifdef IGNORE_RSPACE
                        rsum.clear();
#endif
#ifdef IGNORE_KSPACE
                        ksum.clear();
#endif
                        for (int lm = 0; lm < LEN2; lm++)
                        {
                            gf_r(k, j, i, lm) = rsum.buf[lm] + ksum.buf[lm];
                        }
#if 0
                if(0 == (i|j|k)){
                    rsum.show();
                    ksum.show();
                };
#endif
                    } // for(k, j, i)
        }

        void doFFTOpenImpl()
        {
            MultidimensionalArray<real_t, 3> rbuf;
            MultidimensionalArray<cplx_t, 3> kbuf;
            rbuf.initialize(2 * NZ, 2 * NY, 2 * NX);
            kbuf.initialize(2 * NZ, 2 * NY, 1 + NX);
            fftw_plan plan_fwd =
                fftw_plan_dft_r2c_3d(
                    2 * NZ, 2 * NY, 2 * NX,
                    (double *)(rbuf.getPointer()),
                    (fftw_complex *)(kbuf.getPointer()),
                    FFTW_ESTIMATE);
            for (int lm = 0; lm < LEN2; lm++)
            {
                int i, j, k;
                for (k = 0; k < 2 * NZ; k++)
                    for (int j = 0; j < 2 * NY; j++)
                        for (i = 0; i < 2 * NX; i++)
                        {
                            rbuf(k, j, i) = gf_r(k, j, i, lm);
                        }
                // CALL FFTW
                fftw_execute(plan_fwd);

                for (k = 0; k < 2 * NZ; k++)
                    for (j = 0; j < 2 * NY; j++)
                        for (i = 0; i < 1 + NX; i++)
                        {
                            gf_k(k, j, i, lm) = kbuf(k, j, i);
                        }
            }
            fftw_destroy_plan(plan_fwd);
        }

        void doFFTPeriodicXYZImpl()
        {
            MultidimensionalArray<real_t, 3> rbuf;
            MultidimensionalArray<cplx_t, 3> kbuf;
            rbuf.initialize(NZ, NY, NX);
            kbuf.initialize(NZ, NY, 1 + NX / 2);
            fftw_plan plan_fwd =
                fftw_plan_dft_r2c_3d(
                    NZ, NY, NX,
                    (double *)(rbuf.getPointer()),
                    (fftw_complex *)(kbuf.getPointer()),
                    FFTW_ESTIMATE);
            for (int lm = 0; lm < LEN2; lm++)
            {
                int i, j, k;
                for (k = 0; k < NZ; k++)
                    for (int j = 0; j < NY; j++)
                        for (i = 0; i < NX; i++)
                        {
                            rbuf(k, j, i) = gf_r(k, j, i, lm);
                        }
                // CALL FFTW
                fftw_execute(plan_fwd);

                for (k = 0; k < NZ; k++)
                    for (j = 0; j < NY; j++)
                        for (i = 0; i < 1 + NX / 2; i++)
                        {
                            gf_k(k, j, i, lm) = kbuf(k, j, i);
                        }
            }
            fftw_destroy_plan(plan_fwd);
        }

    public:
        void init(const int bc,
                  const int NX,
                  const int NY,
                  const int NZ)
        {
            this->bc = bc;
            this->NX = NX;
            this->NY = NY;
            this->NZ = NZ;
            if (bc == BOUNDARY_CONDITION_OPEN)
            {
                initOpenImpl();
            }
            else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ)
            {
                initPeriodicXYZImpl();
            }
        }

        void set(const int icut,
                 const double cell_length,
                 const double alpha = 2.4,
                 const int NMAX = 3,
                 const int MMAX = 5)
        {
            if (bc == BOUNDARY_CONDITION_OPEN)
            {
                setOpenImpl(icut, cell_length);
            }
            else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ)
            {
                setPeriodicXYZImpl(icut, cell_length, alpha, NMAX, MMAX);
            }
        }

        void doFFT()
        {
            if (bc == BOUNDARY_CONDITION_OPEN)
            {
                doFFTOpenImpl();
            }
            else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ)
            {
                doFFTPeriodicXYZImpl();
            }
        }

        void transform(const MultidimensionalArray<cplx_t, 4> &mm_k,
                       MultidimensionalArray<cplx_t, 4> &le_k) const
        {
            const int kmax = gf_k.getSize(0);
            const int jmax = gf_k.getSize(1);
            const int imax = gf_k.getSize(2);
            for (int k = 0; k < kmax; k++)
            {
                for (int j = 0; j < jmax; j++)
                {
                    for (int i = 0; i < imax; i++)
                    {
                        typedef Slm<2 * p, cplx_t> slm_t;
                        typedef MultipoleMoment<p, cplx_t> mm_t;
                        typedef LocalExpansion<p, cplx_t> le_t;
                        ((slm_t *)(&gf_k(k, j, i, 0)))
                            ->template transform_M2L<p, p, false>((mm_t *)&mm_k(k, j, i, 0),
                                                                  (le_t *)&le_k(k, j, i, 0));
                    }
                }
            }
        }
    };

    template <int p, int dim>
    static void M2LConvolutionOpenImpl(const GreenFunction<p> &gf,
                                       MultidimensionalArray<Cell_FMM<p>, dim> &cell)
    {
        typedef typename GreenFunction<p>::real_t real_t;
        typedef typename GreenFunction<p>::cplx_t cplx_t;
        enum
        {
            LEN = lmbuf<p>::length,
        };
        // Get # of grid points
        const int NX = gf.NX;
        const int NY = gf.NY;
        const int NZ = gf.NZ;
        // Multipole Moments
        MultidimensionalArray<cplx_t, 4> mm_k;
        // Local Expansions
        MultidimensionalArray<cplx_t, 4> le_k;
        // FFT buffer
        MultidimensionalArray<real_t, 3> rbuf;
        MultidimensionalArray<cplx_t, 3> kbuf;
        // Memory allocation
        mm_k.initialize(2 * NZ, 2 * NY, 1 + NX, LEN);
        le_k.initialize(2 * NZ, 2 * NY, 1 + NX, LEN);
        rbuf.initialize(2 * NZ, 2 * NY, 2 * NX);
        kbuf.initialize(2 * NZ, 2 * NY, 1 + NX);
        // Create plans of FFTW
        fftw_plan plan_fwd =
            fftw_plan_dft_r2c_3d(
                2 * NZ, 2 * NY, 2 * NX,
                (double *)(rbuf.getPointer()),
                (fftw_complex *)(kbuf.getPointer()),
                FFTW_ESTIMATE);
        fftw_plan plan_bkw =
            fftw_plan_dft_c2r_3d(
                2 * NZ, 2 * NY, 2 * NX,
                (fftw_complex *)(kbuf.getPointer()),
                (double *)(rbuf.getPointer()),
                FFTW_ESTIMATE);

        int i, j, k;
        // clear rbuf
        for (k = 0; k < 2 * NZ; k++)
            for (j = 0; j < 2 * NY; j++)
                for (i = 0; i < 2 * NX; i++)
                {
                    rbuf(k, j, i) = 0.0;
                }
        // forward multipole
        for (int lm = 0; lm < LEN; lm++)
        {
            for (k = 0; k < NZ; k++)
                for (j = 0; j < NY; j++)
                    for (i = 0; i < NX; i++)
                    {
                        rbuf(k, j, i) = cell(k, j, i).mm.buf[lm];
                    }

            fftw_execute(plan_fwd);

            for (k = 0; k < 2 * NZ; k++)
                for (j = 0; j < 2 * NY; j++)
                    for (i = 0; i < 1 + NX; i++)
                    {
                        mm_k(k, j, i, lm) = kbuf(k, j, i);
                    }
        }
        // M2L transformation
        gf.transform(mm_k, le_k);
        // backward local expansion
        for (int lm = 0; lm < LEN; lm++)
        {
            for (k = 0; k < 2 * NZ; k++)
                for (j = 0; j < 2 * NY; j++)
                    for (i = 0; i < 1 + NX; i++)
                    {
                        kbuf(k, j, i) = le_k(k, j, i, lm);
                    }

            fftw_execute(plan_bkw);

            const double norm = 1.0 / (8 * NX * NY * NZ);
            for (k = 0; k < NZ; k++)
                for (j = 0; j < NY; j++)
                    for (i = 0; i < NX; i++)
                    {
                        cell(k, j, i).le.buf[lm] = norm * rbuf(k, j, i);
                    }
        }
        fftw_destroy_plan(plan_fwd);
        fftw_destroy_plan(plan_bkw);
    }

    template <int p, int dim>
    static void M2LConvolutionPeriodicXYZImpl(const GreenFunction<p> &gf,
                                              MultidimensionalArray<Cell_FMM<p>, dim> &cell)
    {
        typedef typename GreenFunction<p>::real_t real_t;
        typedef typename GreenFunction<p>::cplx_t cplx_t;
        enum
        {
            LEN = lmbuf<p>::length,
        };
        // Get # of grid points
        const int NX = gf.NX;
        const int NY = gf.NY;
        const int NZ = gf.NZ;
        // Multipole Moments
        MultidimensionalArray<cplx_t, 4> mm_k;
        // Local Expansions
        MultidimensionalArray<cplx_t, 4> le_k;
        // FFT buffer
        MultidimensionalArray<real_t, 3> rbuf;
        MultidimensionalArray<cplx_t, 3> kbuf;
        // Memory allocation
        mm_k.initialize(NZ, NY, 1 + NX / 2, LEN);
        le_k.initialize(NZ, NY, 1 + NX / 2, LEN);
        rbuf.initialize(NZ, NY, NX);
        kbuf.initialize(NZ, NY, 1 + NX / 2);
        // Create plans of FFTW
        fftw_plan plan_fwd =
            fftw_plan_dft_r2c_3d(
                NZ, NY, NX,
                (double *)(rbuf.getPointer()),
                (fftw_complex *)(kbuf.getPointer()),
                FFTW_ESTIMATE);
        fftw_plan plan_bkw =
            fftw_plan_dft_c2r_3d(
                NZ, NY, NX,
                (fftw_complex *)(kbuf.getPointer()),
                (double *)(rbuf.getPointer()),
                FFTW_ESTIMATE);

        int i, j, k;
        // forward multipole
        for (int lm = 0; lm < LEN; lm++)
        {
            for (k = 0; k < NZ; k++)
                for (j = 0; j < NY; j++)
                    for (i = 0; i < NX; i++)
                    {
                        rbuf(k, j, i) = cell(k, j, i).mm.buf[lm];
                    }

            fftw_execute(plan_fwd);

            for (k = 0; k < NZ; k++)
                for (j = 0; j < NY; j++)
                    for (i = 0; i < 1 + NX / 2; i++)
                    {
                        mm_k(k, j, i, lm) = kbuf(k, j, i);
                    }
        }
        // M2L transformation
        gf.transform(mm_k, le_k);
        // backward local expansion
        for (int lm = 0; lm < LEN; lm++)
        {
            for (k = 0; k < NZ; k++)
                for (j = 0; j < NY; j++)
                    for (i = 0; i < 1 + NX / 2; i++)
                    {
                        kbuf(k, j, i) = le_k(k, j, i, lm);
                    }

            fftw_execute(plan_bkw);

            const double norm = 1.0 / (NX * NY * NZ);
            for (k = 0; k < NZ; k++)
                for (j = 0; j < NY; j++)
                    for (i = 0; i < NX; i++)
                    {
                        cell(k, j, i).le.buf[lm] = norm * rbuf(k, j, i);
                    }
        }
        fftw_destroy_plan(plan_fwd);
        fftw_destroy_plan(plan_bkw);
    }

    template <int p, int dim>
    void M2LConvolution(const GreenFunction<p> &gf,
                        MultidimensionalArray<Cell_FMM<p>, dim> &cell,
                        const int bc)
    {
        if (bc == BOUNDARY_CONDITION_OPEN)
        {
            M2LConvolutionOpenImpl<p, dim>(gf, cell);
        }
        else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ)
        {
            M2LConvolutionPeriodicXYZImpl<p, dim>(gf, cell);
        }
    }

} // END of namespace of ParticleSimulator
