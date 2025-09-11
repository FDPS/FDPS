/*
 Derived from PMMM (Apache License 2.0)
 Original Copyright (c) 2012-2014 Keigo Nitadori
 Modifications (c) 2025 FDPS team
 Licensed under the Apache License, Version 2.0 (see LICENSE)
*/
#pragma once
#include <iostream>
#include <complex>
#include <algorithm>
#include <cassert>
#include <unordered_map>

namespace ParticleSimulator {

inline double factorial(const int i) {
    assert(i >= 0);
    return i ? double(i) * factorial(i - 1) : 1.0;
}

inline double factinv(const int i) { return 1.0 / factorial(i); }

template <class T>
inline T tp_zero_val() {
    return T(0);
}

template <class real_t = double, class cplx_t = std::complex<real_t> >
class lmbuf {
   private:
    enum class status_t : int { is_default = 0, is_reinterpreting = 1 };

   public:
    status_t stat;
    int p;
    real_t *buf;

    lmbuf() : stat(status_t::is_default), p(-1), buf(nullptr) {}

    ~lmbuf() { this->freeMem(); }

    lmbuf(const lmbuf &src) {  // copy constructor
        this->stat = status_t::is_default;
        // [Note (tag: #55334c30)]
        //    Memory used by a copy of lmbuf class must be released
        //    after use. Hence, we should set status_t::is_default
        //    this->stat regardless of the value of src.stat.
        this->p = src.p;
        if (src.buf != nullptr) {
            assert(src.p >= 0);
            const int length = (src.p + 1) * (src.p + 1);
            buf = new real_t[length];
            for (int i = 0; i < length; i++) this->buf[i] = src.buf[i];
        } else {
            this->buf = nullptr;
        }
    }

    void alloc(int p, const bool clear_flag = false) {
        // This function is used to allocate memory for buf
        // and can be used if and only if this->stat ==
        // status_t::is_default && p >= 0.
        assert(this->stat == status_t::is_default);
        assert(p >= 0);
        if (this->p != p) {
            // In this case, buf == nullptr or the size of buf
            // is not consistent with that required by the specified p.
            // So, we (re)allocate memory for buf.
            this->p = p;
            const int length = (p + 1) * (p + 1);
            if (buf != nullptr) delete[] buf;
            // 20250815 work for monopole moment
            // buf = new real_t[length];
            if (p == 0) {
                buf = new real_t[4];  // even if monopole moment, dipole correction is required.
            } else {
                buf = new real_t[length];
            }

        } else {
            // In this case, we reuse buf.
            assert(buf != nullptr);
        }
        if (clear_flag) this->clear();
    }

    void reinterpret(int p, real_t *adr_buf) {
        // This function is used to use an existing array
        // as an instance of lmbuf class;
        if (this->stat == status_t::is_default) {  // to check memory leak.
            assert(buf == nullptr);
        }
        this->stat = status_t::is_reinterpreting;
        this->p = p;
        buf = adr_buf;
        static real_t prefetch = buf[0];
    }

    void moveFrom(lmbuf<real_t, cplx_t> &src) {
        assert(p == src.p);
        buf = src.buf;
        src.buf = nullptr;
    }

    void freeMem() {
        // This function explicitly frees memory if
        // buf == status_t::is_default.
        if (stat == status_t::is_default) {
            if (buf != nullptr) {
                delete[] buf;
                buf = nullptr;
                p = -1;
            }
        }
    }

    void clear() {
        if (buf != nullptr) {
            const int length = (p + 1) * (p + 1);
            for (int n = 0; n < length; n++) buf[n] = tp_zero_val<real_t>();
        }
    }

    int size() const { return (p + 1) * (p + 1); }

    // X_l^{-m} = (-1)^m conj(X_l^m)
    cplx_t val_at(const int l, const int m) const __restrict {
#ifdef PARTICLE_SIMULATOR_RANGE_CHECK_LMBUF
        assert(0 <= l && l <= p);
        assert(abs(m) <= l);
#endif
        const real_t *base = &buf[0] + l * (l + 1);

        real_t re = base[+abs(m)];
        real_t im = base[-abs(m)];

        if (m == 0) im = tp_zero_val<real_t>();
        if (m < 0) {
            if (m & 1)
                re = -re;
            else
                im = -im;
        }
        return cplx_t(re, im);
    }
    void put_at(const int l, const int m, const cplx_t z) {
#ifdef PARTICLE_SIMULATOR_RANGE_CHECK_LMBUF
        assert(0 <= l && l <= p);
        assert(0 <= m && m <= l);
#endif
        /*     */ buf[l * (l + 1) + m] = real(z);
        if (m > 0) buf[l * (l + 1) - m] = imag(z);
    }
    void put_at(const int l, const int m, const real_t re, const real_t im) {
#ifdef PARTICLE_SIMULATOR_RANGE_CHECK_LMBUF
        assert(0 <= l && l <= p);
        assert(0 <= m && m <= l);
#endif
        /*     */ buf[l * (l + 1) + m] = re;
        if (m > 0) buf[l * (l + 1) - m] = im;
    }
    void accum_at(const int l, const int m, const cplx_t z) {
#ifdef PARTICLE_SIMULATOR_RANGE_CHECK_LMBUF
        assert(0 <= l && l <= p);
        assert(0 <= m && m <= l);
#endif
        /*     */ buf[l * (l + 1) + m] += real(z);
        if (m > 0) buf[l * (l + 1) - m] += imag(z);
    }
    void accum_at(const int l, const int m, const real_t re, const real_t im) {
#ifdef PARTICLE_SIMULATOR_RANGE_CHECK_LMBUF
        assert(0 <= l && l <= p);
        assert(0 <= m && m <= l);
#endif
        /*     */ buf[l * (l + 1) + m] += re;
        if (m > 0) buf[l * (l + 1) - m] += im;
    }

    // only for positive m
    static int idx_re(const int l, const int m) { return l * (l + 1) + m; }
    static int idx_im(const int l, const int m) { return l * (l + 1) - m; }
    real_t get_re(const int l, const int m) const { return buf[idx_re(l, m)]; }
    real_t get_im(const int l, const int m) const { return buf[idx_im(l, m)]; }

    void show(FILE *fp = stdout, const char *fmt = "%+f %+fi, ") const {
        for (int l = 0; l <= p; l++) {
            fprintf(fp, "l=%d: ", l);
            for (int m = 0; m <= l; m++) {
                const cplx_t z = val_at(l, m);
                fprintf(fp, fmt, real(z), imag(z));
            }
            fprintf(fp, "\n");
        }
    }
    void show_diff(const lmbuf &rhs, FILE *fp = stdout, const char *fmt = "%+e %+ei, ") const {
        for (int l = 0; l <= p; l++) {
            fprintf(fp, "l=%d: ", l);
            for (int m = 0; m <= l; m++) {
#if 1
                const cplx_t z = val_at(l, m) - rhs.val_at(l, m);
#else
                cplx_t z = val_at(l, m) - rhs.val_at(l, m);
                z /= abs(val_at(l, m));
#endif
                fprintf(fp, fmt, real(z), imag(z));
            }
            fprintf(fp, "\n");
        }
    }
    void show_all(FILE *fp = stdout, const char *fmt = "%+f %+fi, ") const {
        for (int l = 0; l <= p; l++) {
            fprintf(fp, "l=%d: ", l);
            for (int m = -l; m <= l; m++) {
                const cplx_t z = val_at(l, m);
                fprintf(fp, fmt, real(z), imag(z));
            }
            fprintf(fp, "\n");
        }
    }
};

// Rlm = (r^l / (l+m)!) P_l^m (cos\theta) exp(im\phi), for m>=0
template <class real_t = double, class cplx_t = std::complex<real_t> >
struct Rlm : public lmbuf<real_t, cplx_t> {
    static int p_limit;
    static std::vector<real_t> tbl_inv;
    static std::vector<real_t> tbl_factinv;

    void make_table(const int p_in) {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp critical
#endif
        {
            if (p_in > p_limit) {
                p_limit = p_in;
                tbl_inv.resize(p_limit + 1);
                tbl_inv[0] = real_t(1.0);  // dummy value; tbl_inv[0] is not used.
                for (int i = 1; i < p_limit + 1; i++) {
                    tbl_inv[i] = real_t(1.0) / real_t(i);
                }
                tbl_factinv.resize(2 * p_limit + 1);
                for (int i = 0; i < 2 * p_limit + 1; i++) {
                    assert(factorial(i) > 0);
                    tbl_factinv[i] = 1.0 / factorial(i);
                }
            }
        }
    }

    void eval_opt(const real_t x, const real_t y, const real_t z) {
        assert(p_limit >= this->p);

        const real_t r2 = x * x + y * y + z * z;
        const cplx_t u(x, y);
        const real_t zz = z;

        cplx_t um(1.0, 0.0);
        real_t pmm = 1.0;
        real_t _2mp1 = 1.0;
        for (int m = 0; m <= this->p; m++) {
            real_t p0 = pmm;
            real_t p1 = _2mp1 * zz * pmm;
            real_t plm;
            for (int l = m; l <= this->p;) {
                plm = p0;
                const cplx_t val = (tbl_factinv[l + m] * plm) * um;
                this->put_at(l, m, val);
                break;
            }
            for (int l = m + 1; l <= this->p;) {
                plm = p1;
                const cplx_t val = (tbl_factinv[l + m] * plm) * um;
                this->put_at(l, m, val);
                break;
            }
            real_t c0 = _2mp1;
            real_t c1 = c0 + 2.0;
            for (int l = m + 2; l <= this->p; l++) {
                plm = (tbl_inv[l - m]) * (c1 * zz * p1 - c0 * r2 * p0);
                p0 = p1;
                p1 = plm;
                c0 += 1.0;
                c1 += 2.0;
                const cplx_t val = (tbl_factinv[l + m] * plm) * um;
                this->put_at(l, m, val);
            }
            // end of for(m), update values
            um *= u;
            pmm *= real_t(-(2 * m + 1));
            _2mp1 += 2.0;
        }
    }
};
template <class real_t, class cplx_t>
int Rlm<real_t, cplx_t>::p_limit = -1;
template <class real_t, class cplx_t>
std::vector<real_t> Rlm<real_t, cplx_t>::tbl_inv;
template <class real_t, class cplx_t>
std::vector<real_t> Rlm<real_t, cplx_t>::tbl_factinv;

// Slm = (-1)^{l+m} ((l-m)! / r^{l+1}) P_l^m (cos\theta) exp(im\phi), for m>=0
template <class real_t = double, class cplx_t = std::complex<real_t> >
struct Slm : public lmbuf<real_t, cplx_t> {
    static int p_limit;
    static std::vector<real_t> tbl_inv;

    void make_table(const int p_in) {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp critical
#endif
        {
            if (p_in > p_limit) {
                p_limit = p_in;
                tbl_inv.resize(p_limit + 1);
                tbl_inv[0] = real_t(1.0);  // dummy value; tbl_inv[0] is not used.
                for (int i = 1; i < p_limit + 1; i++) {
                    tbl_inv[i] = real_t(1.0) / real_t(i);
                }
            }
        }
    }

    void eval_opt(const real_t x, const real_t y, const real_t z) {
        assert(p_limit >= this->p);

        const real_t r2 = x * x + y * y + z * z;
        const cplx_t u(x, y);  // for (-1)^(l+m)
        const real_t zz = -z;  // for (-1)^l

        const real_t rinv2 = real_t(1.0) / r2;
        const real_t rinv = sqrt(rinv2);

        cplx_t um(1.0, 0.0);
        real_t pmm = 1.0;
        real_t rinv_2mp1 = rinv;
        real_t _2mp1 = 1.0;
        for (int m = 0; m <= this->p; m++) {
            real_t rinv_2lp1 = rinv_2mp1;

            real_t p0, p1, plm;
            p0 = p1 = 0.0;  // avoid warning
            for (int l = m; l <= this->p;) {
                plm = p0 = pmm;
                const cplx_t val = (rinv_2lp1 * plm) * um;
                this->put_at(l, m, val);
                rinv_2lp1 *= rinv2;
                break;
            }
            for (int l = m + 1; l <= this->p;) {
                plm = p1 = _2mp1 * zz * pmm;
                const cplx_t val = (rinv_2lp1 * plm) * um;
                this->put_at(l, m, val);
                rinv_2lp1 *= rinv2;
                break;
            }
            real_t fact_lmm = 2.0;
            real_t c0 = _2mp1;
            real_t c1 = c0 + 2.0;
            for (int l = m + 2; l <= this->p; l++) {
                plm = (tbl_inv[l - m]) * (c1 * zz * p1 - c0 * r2 * p0);
                p0 = p1;
                p1 = plm;
                c0 += 1.0;
                c1 += 2.0;
                const cplx_t val = (rinv_2lp1 * fact_lmm * plm) * um;
                this->put_at(l, m, val);
                rinv_2lp1 *= rinv2;
                fact_lmm *= real_t(l - m + 1);
            }
            // end of for(m), update values
            um *= u;
            pmm *= -_2mp1;
            rinv_2mp1 *= rinv2;
            _2mp1 += 2.0;
        }
    }

    void transform_M2L(const lmbuf<real_t> &mm, lmbuf<real_t> &le, const bool accumulate) const __restrict {
        const int psrc = mm.p;
        const int pdst = le.p;
        assert(this->p >= psrc + pdst);
        for (int ldst = 0; ldst <= pdst; ldst++) {
            for (int mdst = 0; mdst <= ldst; mdst++) {
                real_t rsum = 0.0;
                real_t isum = 0.0;
                for (int lsrc = 0; lsrc <= psrc; lsrc++) {
                    do {
                        const int msrc = 0;
                        const real_t mr = mm.get_re(lsrc, msrc);
                        const cplx_t gp = this->val_at(lsrc + ldst, msrc + mdst);
                        const real_t a = real(gp);
                        const real_t b = imag(gp);
                        rsum += a * mr;
                        isum += b * mr;
                    } while (0);
                    for (int msrc = 1; msrc <= lsrc; /* msrc++ */) {
                        {
                            const real_t mr = mm.get_re(lsrc, msrc);
                            const real_t mi = mm.get_im(lsrc, msrc);
                            const cplx_t gp = this->val_at(lsrc + ldst, mdst + msrc);
                            const cplx_t gm = this->val_at(lsrc + ldst, mdst - msrc);
                            const real_t a = real(gp);
                            const real_t b = imag(gp);
                            const real_t c = real(gm);
                            const real_t d = imag(gm);
                            rsum += (a - c) * mr - (b + d) * mi;
                            isum += (b - d) * mr + (a + c) * mi;
                        }
                        if (++msrc > lsrc) break;
                        {
                            const real_t mr = mm.get_re(lsrc, msrc);
                            const real_t mi = mm.get_im(lsrc, msrc);
                            const cplx_t gp = this->val_at(lsrc + ldst, mdst + msrc);
                            const cplx_t gm = this->val_at(lsrc + ldst, mdst - msrc);
                            const real_t a = real(gp);
                            const real_t b = imag(gp);
                            const real_t c = real(gm);
                            const real_t d = imag(gm);
                            rsum += (a + c) * mr - (b - d) * mi;
                            isum += (b + d) * mr + (a - c) * mi;
                        }
                        if (++msrc > lsrc) break;
                    }  // for(msrc)
                }  // for(lsrc)
                if (accumulate) {
                    le.accum_at(ldst, mdst, cplx_t(rsum, isum));
                } else {
                    le.put_at(ldst, mdst, cplx_t(rsum, isum));
                }
            }
        }
    }
};
template <class real_t, class cplx_t>
int Slm<real_t, cplx_t>::p_limit = -1;
template <class real_t, class cplx_t>
std::vector<real_t> Slm<real_t, cplx_t>::tbl_inv;

template <class real_t = double, class cplx_t = std::complex<real_t> >
struct MultipoleMoment : public lmbuf<real_t, cplx_t> {
    // P2M
    void assign_particle(const real_t x, const real_t y, const real_t z, const real_t charge) {
        Rlm<real_t, cplx_t> rlm;
        rlm.alloc(this->p);
        rlm.eval_opt(x, y, z);  // generate R_l^{m}
        const int len = (this->p + 1) * (this->p + 1);
        for (int n = 0; n < len; n++) {
            this->buf[n] += charge * rlm.buf[n];
        }
        rlm.freeMem();
    }
    void assign_particle(const Vector3<real_t> &pos_this, const Vector3<real_t> &pos_ptcl, const real_t charge) {
        const Vector3<real_t> dr = pos_this - pos_ptcl;
        assign_particle(dr.x, dr.y, dr.z, charge);
    }
    // M2P
    real_t eval_potential(const real_t x, const real_t y, const real_t z) const {
        Slm<real_t, cplx_t> slm;
        slm.alloc(this->p);
        slm.eval_opt(-x, y, z);  // generate S_l^{-m}

        real_t pot = 0.0;
        for (int l = 0; l <= this->p; l++) {
            pot += real(this->val_at(l, 0) * slm.val_at(l, 0));
            for (int m = 1; m <= l; m++) {
                pot += 2.0 * real(this->val_at(l, m) * slm.val_at(l, m));
            }
        }
        slm.freeMem();
        return pot;
    }
    real_t eval_potential(const Vector3<real_t> &pos_this, const Vector3<real_t> &pos_ptcl) const {
        const Vector3<real_t> dr = pos_ptcl - pos_this;
        return eval_potential(dr.x, dr.y, dr.z);
    }

    // M2M
    void assign_from_MM(const MultipoleMoment<real_t, cplx_t> &MMsrc, const real_t x, const real_t y, const real_t z) {
        const int psrc = MMsrc.p;
        Rlm<real_t, cplx_t> rlm;
        rlm.alloc(psrc);
        rlm.eval_opt(x, y, z);  // R_{l}^{m}

        for (int l = 0; l <= this->p; l++) {
            for (int m = 0; m <= l; m++) {
                cplx_t val(0.0, 0.0);
                // for(int lsrc=0; lsrc<=l; lsrc++){
                for (int lsrc = 0; lsrc <= std::min(l, psrc); lsrc++) {
                    const int mbeg = std::max(-lsrc, m - (l - lsrc));
                    const int mend = std::min(+lsrc, m + (l - lsrc));
                    for (int msrc = mbeg; msrc <= mend; msrc++) {
                        val += MMsrc.val_at(l - lsrc, m - msrc) * rlm.val_at(lsrc, msrc);
                    }
                }
                this->accum_at(l, m, val);
            }
        }
        rlm.freeMem();
    }
    void assign_from_MM(const MultipoleMoment<real_t, cplx_t> &MMsrc, const Vector3<real_t> &pos_this, const Vector3<real_t> &pos_src) {
        const Vector3<real_t> dr = pos_this - pos_src;
        this->assign_from_MM(MMsrc, dr.x, dr.y, dr.z);
    }
};

template <class real_t = double, class cplx_t = std::complex<real_t> >
struct LocalExpansion : public lmbuf<real_t, cplx_t> {
    void read_phi_and_grad(real_t &phi, real_t &ex, real_t &ey, real_t &ez) const {
        phi = +this->buf[0];  // (0,0)
        ex = -this->buf[3];   // Re(1,1)
        ey = +this->buf[1];   // Im(1,1)
        ez = +this->buf[2];   // (1,0)
    }
    void read_phi_and_grad(real_t &phi, Vector3<real_t> &v) const { read_phi_and_grad(phi, v.x, v.y, v.z); }

    // L2P
    real_t eval_potential(const real_t x, const real_t y, const real_t z) const {
        Rlm<real_t, cplx_t> rlm;
        rlm.alloc(this->p);
        rlm.eval_opt(x, y, z);

        real_t pot = 0.0;
        for (int l = 0; l <= this->p; l++) {
            pot += real(this->val_at(l, 0) * rlm.val_at(l, 0));
            for (int m = 1; m <= l; m++) {
                pot += 2.0 * real(this->val_at(l, m) * rlm.val_at(l, m));
            }
        }
        rlm.freeMem();
        return pot;
    }
    real_t eval_potential(const Vector3<real_t> &pos_this, const Vector3<real_t> &pos_ptcl) const {
        const Vector3<real_t> dr = pos_ptcl - pos_this;
        return eval_potential(dr.x, dr.y, dr.z);
    }

    // L2L
    void assign_from_LE(const LocalExpansion<real_t, cplx_t> &LEsrc, const real_t x, const real_t y, const real_t z) {
        const int psrc = LEsrc.p;
        Rlm<real_t, cplx_t> rlm;
        rlm.alloc(psrc);
        rlm.eval_opt(x, y, z);

        for (int l = 0; l <= this->p; l++) {
            for (int m = 0; m <= l; m++) {
                cplx_t val(0.0, 0.0);
                for (int lsrc = l; lsrc <= psrc; lsrc++) {
                    for (int msrc = -(lsrc - l) + m; msrc <= (lsrc - l) + m; msrc++) {
                        val += LEsrc.val_at(lsrc, msrc) * rlm.val_at(lsrc - l, msrc - m);
                        /*
                            printf("---\n");
                            printf("L'(%d, %d) <= L(%d, %d) *  R(%d, %d)\n",
                                   l, m, lsrc, msrc, lsrc-l, msrc-m);
                            printf("L = (%e, %e)\n",
                                   LEsrc.val_at(lsrc, msrc).real(),
                                   LEsrc.val_at(lsrc, msrc).imag());
                            printf("R = (%e, %e)\n",
                                   rlm.val_at(lsrc-l, msrc-m).real(),
                                   rlm.val_at(lsrc-l, msrc-m).imag());
                            printf("val = (%e, %e)\n", val.real(), val.imag());
                        */
                    }
                }
                this->accum_at(l, m, val);
            }
        }
        rlm.freeMem();
    }
    void assign_from_LE(const LocalExpansion<real_t, cplx_t> &LEsrc, const Vector3<real_t> &pos_this, const Vector3<real_t> &pos_src) {
        const Vector3<real_t> dr = pos_this - pos_src;
        this->assign_from_LE(LEsrc, dr.x, dr.y, dr.z);
    }
    void assign_from_LE_opt(const LocalExpansion<real_t, cplx_t> &LEsrc, const real_t x, const real_t y, const real_t z, Rlm<real_t, cplx_t> &rlm) {
        const int psrc = LEsrc.p;
        assert(rlm.p == psrc);
        rlm.eval_opt(x, y, z);

        for (int l = 0; l <= this->p; l++) {
            for (int m = 0; m <= l; m++) {
                cplx_t val(0.0, 0.0);
                for (int lsrc = l; lsrc <= psrc; lsrc++) {
                    for (int msrc = -(lsrc - l) + m; msrc <= (lsrc - l) + m; msrc++) {
                        val += LEsrc.val_at(lsrc, msrc) * rlm.val_at(lsrc - l, msrc - m);
                        /*
                            printf("---\n");
                            printf("L'(%d, %d) <= L(%d, %d) *  R(%d, %d)\n",
                                   l, m, lsrc, msrc, lsrc-l, msrc-m);
                            printf("L = (%e, %e)\n",
                                   LEsrc.val_at(lsrc, msrc).real(),
                                   LEsrc.val_at(lsrc, msrc).imag());
                            printf("R = (%e, %e)\n",
                                   rlm.val_at(lsrc-l, msrc-m).real(),
                                   rlm.val_at(lsrc-l, msrc-m).imag());
                            printf("val = (%e, %e)\n", val.real(), val.imag());
                        */
                    }
                }
                this->accum_at(l, m, val);
            }
        }
    }
    void assign_from_LE_opt(const LocalExpansion<real_t, cplx_t> &LEsrc, const Vector3<real_t> &pos_this, const Vector3<real_t> &pos_src,
                            Rlm<real_t, cplx_t> &rlm) {
        const Vector3<real_t> dr = pos_this - pos_src;
        this->assign_from_LE_opt(LEsrc, dr.x, dr.y, dr.z, rlm);
    }
    // M2L
    void assign_from_MM(const MultipoleMoment<real_t, cplx_t> &MMsrc, const real_t x, const real_t y, const real_t z) {
        const int psrc = MMsrc.p;
        Slm<real_t, cplx_t> slm;
        slm.alloc(this->p + psrc);
        slm.eval_opt(-x, y, z);  // eval S_l^{-m}

        for (int l = 0; l <= this->p; l++) {
            for (int m = 0; m <= l; m++) {
                cplx_t val(0.0, 0.0);
                for (int lsrc = 0; lsrc <= psrc; lsrc++) {
                    for (int msrc = -lsrc; msrc <= lsrc; msrc++) {
                        val += MMsrc.val_at(lsrc, msrc) * slm.val_at(lsrc + l, msrc + m);
                        // printf("L(%d, %d) <= M(%d, %d) *  L(%d, %d)\n",
                        //          l, m, lsrc, msrc, lsrc+l, msrc-m);
                    }
                }
                this->accum_at(l, m, val);
            }
        }
        slm.freeMem();
    }
    void assign_from_MM(const MultipoleMoment<real_t, cplx_t> &MMsrc, const Vector3<real_t> &pos_this, const Vector3<real_t> &pos_src) {
        const Vector3<real_t> dr = pos_this - pos_src;
        this->assign_from_MM(MMsrc, dr.x, dr.y, dr.z);
    }

    void assign_from_MM_opt(const MultipoleMoment<real_t, cplx_t> &MMsrc, const real_t x, const real_t y, const real_t z) {
        const int psrc = MMsrc.p;
        Slm<real_t, cplx_t> slm;
        slm.alloc(this->p + psrc);
        slm.eval_opt(-x, y, z);  // eval S_l^{-m}

        for (int ldst = 0; ldst <= this->p; ldst++) {
            for (int mdst = 0; mdst <= ldst; mdst++) {
                real_t rsum = 0.0;
                real_t isum = 0.0;
                for (int lsrc = 0; lsrc <= psrc; lsrc++) {
                    do {
                        const int msrc = 0;
                        const real_t mr = MMsrc.get_re(lsrc, msrc);
                        const cplx_t gp = slm.val_at(lsrc + ldst, msrc + mdst);
                        const real_t a = real(gp);
                        const real_t b = imag(gp);
                        rsum += a * mr;
                        isum += b * mr;
                    } while (0);
                    for (int msrc = 1; msrc <= lsrc; /* msrc++ */) {
                        {
                            const real_t mr = MMsrc.get_re(lsrc, msrc);
                            const real_t mi = MMsrc.get_im(lsrc, msrc);
                            const cplx_t gp = slm.val_at(lsrc + ldst, mdst + msrc);
                            const cplx_t gm = slm.val_at(lsrc + ldst, mdst - msrc);
                            const real_t a = real(gp);
                            const real_t b = imag(gp);
                            const real_t c = real(gm);
                            const real_t d = imag(gm);
                            rsum += (a - c) * mr - (b + d) * mi;
                            isum += (b - d) * mr + (a + c) * mi;
                        }
                        if (++msrc > lsrc) break;
                        {
                            const real_t mr = MMsrc.get_re(lsrc, msrc);
                            const real_t mi = MMsrc.get_im(lsrc, msrc);
                            const cplx_t gp = slm.val_at(lsrc + ldst, mdst + msrc);
                            const cplx_t gm = slm.val_at(lsrc + ldst, mdst - msrc);
                            const real_t a = real(gp);
                            const real_t b = imag(gp);
                            const real_t c = real(gm);
                            const real_t d = imag(gm);
                            rsum += (a + c) * mr - (b - d) * mi;
                            isum += (b + d) * mr + (a - c) * mi;
                        }
                        if (++msrc > lsrc) break;
                    }  // for(msrc)
                }  // for(lsrc)
                this->accum_at(ldst, mdst, cplx_t(rsum, isum));
            }
        }
        slm.freeMem();
    }
    void assign_from_MM_opt(const MultipoleMoment<real_t, cplx_t> &MMsrc, const Vector3<real_t> &pos_this, const Vector3<real_t> &pos_src) {
        const Vector3<real_t> dr = pos_this - pos_src;
        this->assign_from_MM_opt(MMsrc, dr.x, dr.y, dr.z);
    }
};

}  // namespace ParticleSimulator
