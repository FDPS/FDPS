#pragma once
struct Force {
    PS::F64vec grad;
    PS::F64 pot;
    void clear() {
        grad = 0.0;
        pot = 0.0;
    }
    void addForcePMMM(const PS::F64vec &f, const PS::F64 p) {
        grad += f;
        pot += p;
    }
};

struct FP {
    PS::S64 id;
    PS::F64 charge;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec grad;
    PS::F64 pot;
    PS::F64vec getPos() const { return pos; }
    PS::F64 getCharge() const { return charge; }
    PS::F64 getChargeParticleMesh() const { return charge; }
    void copyFromForce(const Force &f) {
        this->grad = f.grad;
        this->pot = f.pot;
    }
    void copyFromForceParticleMesh(const PS::F64vec &g) { this->grad = g; }
    void copyFromForcePMMM(const Force &f) {  // copy force from PMMM
        this->grad = f.grad;
        this->pot = f.pot;
    }
    void setPos(const PS::F64vec &pos_new) { this->pos = pos_new; }
    void clear() {
        this->grad = 0.0;
        this->pot = 0.0;
    }
    void dump(std::ostream &fout = std::cout) {
        std::cout << "id= " << id << " charge= " << charge << " pos= " << pos << " vel= " << vel << " grad= " << grad << " pot= " << pot << std::endl;
    }
};

struct EP {
    PS::F64 charge;
    PS::F64vec pos;
    static inline PS::F64 r_search;
    PS::F64vec getPos() const { return pos; }
    void setPos(const PS::F64vec &pos_new) { this->pos = pos_new; }
    PS::F64 getCharge() const { return charge; }
    PS::F64 getRSearch() const { return r_search; }
    void copyFromFP(const FP &fp) {
        this->charge = fp.charge;
        this->pos = fp.pos;
    }
};

struct EwaldEngine {
    static constexpr PS::F64 pi = 3.14159265358979323846;
    static constexpr PS::F64 inv_sqrt_pi = 1.0 / sqrt(pi);
    static inline PS::F64 q_tot;
    static inline PS::F64 alpha;
    static inline PS::F64vec len_unit_cell;
    static inline PS::F64 volume_unit_cell;
    static inline PS::S32 n_img_r = 3;  // Ewald image number in real space
    static inline PS::S32 n_img_k = 5;  // Ewald image number in reciprocal space
    static inline PS::F64vec ***RLV;    // Reciprocal Lattice Vector
    static inline PS::F64 ***qcos_glb;  // qcos = sum (q * cos(RLV * pos))
    static inline PS::F64 ***qsin_glb;  // qsin = sum (q * sin(RLV * pos))
    static inline PS::F64vec *RLV_buf;
    static inline PS::F64 *qcos_glb_buf;
    static inline PS::F64 *qsin_glb_buf;

    template <typename Tpsys>
    static void setQTot(Tpsys &ptcl) {
        PS::F64 q_tot_loc = 0.0;
        for (PS::S32 i = 0; i < ptcl.getNumberOfParticleLocal(); i++) {
            q_tot_loc += ptcl[i].charge;
        }
        q_tot = PS::Comm::getSum(q_tot_loc);
    }
    static void setAlpha(PS::F64 alpha_new) { alpha = alpha_new; }
    static void setUnitCell(PS::F64ort pos_unit_cell) {
        len_unit_cell = pos_unit_cell.getFullLength();
        volume_unit_cell = len_unit_cell.x * len_unit_cell.y * len_unit_cell.z;
    }
    static void setNImgR(PS::S32 n) { n_img_r = n; }
    static void setNImgK(PS::S32 n) { n_img_k = n; }

    template <typename Tpi, typename Tpj, typename Tforce>
    static void calcForceSR(EP *epi, const PS::S32 n_epi, Tpj *pj, const PS::S32 n_pj, Tforce *force) {
        if (n_pj <= 0) return;  // It is necessary to avoid to adding the potential offset when n_pj = 0
        for (PS::S32 i = 0; i < n_epi; i++) {
            PS::F64vec grad = 0.0;
            PS::F64 pot = 0.0;
            const PS::F64vec posi = epi[i].pos;
            for (int j = 0; j < n_pj; j++) {
                const PS::F64vec posj = pj[j].getPos();
                const PS::F64vec rij = posi - posj;
                if (rij.x == 0.0 && rij.y == 0.0 && rij.z == 0.0) continue;
                auto dr2 = rij * rij;
                auto inv_dr = 1.0 / sqrt(dr2);
                auto inv_dr2 = inv_dr * inv_dr;
                auto dr = dr2 * inv_dr;
                auto alphadr = alpha * dr;
                auto qj = pj[j].charge;
                grad -= qj * (erfc(alphadr) * inv_dr + 2.0 * alpha * inv_sqrt_pi * exp(-alphadr * alphadr)) * inv_dr2 * rij;
                pot += qj * erfc(alphadr) * inv_dr;
            }
            force[i].grad += grad;
            force[i].pot += pot;
        }
    }

    template <typename Tpsys>
    static void energyCorrection(Tpsys &ptcl) {
        for (int i = 0; i < ptcl.getNumberOfParticleLocal(); i++) {
            // self energy correction
            // the last factor 2 is for 1/2 in the potential energy formula 1/2 * (q0*pot0 + q1*pot1 + ...) in the Ewald summation
            ptcl[i].pot -= ptcl[i].getCharge() * alpha * inv_sqrt_pi * 2.0;
            // charged system correction
            ptcl[i].pot -= pi * q_tot / (alpha * alpha * volume_unit_cell);
        }
    }

    template <typename Tpsys>
    static void calcForceLR(Tpsys &ptcl) {
        const PS::S32 NIMG = n_img_k;
        const PS::S32 n = ptcl.getNumberOfParticleLocal();
        static bool first_call = true;
        if (first_call) {
            auto NIMG = n_img_k;
            const PS::F64 gx = 2.0 * pi / len_unit_cell.x;
            const PS::F64 gy = 2.0 * pi / len_unit_cell.y;
            const PS::F64 gz = 2.0 * pi / len_unit_cell.z;
            const int S = 2 * NIMG + 1;
            PS::F64 qcos_loc[S][S][S];  // qcos = sum (q * cos(RLV * pos))
            PS::F64 qsin_loc[S][S][S];  // qsin = sum (q * sin(RLV * pos))
            RLV_buf = new PS::F64vec[S * S * S];
            qcos_glb_buf = new PS::F64[S * S * S];
            qsin_glb_buf = new PS::F64[S * S * S];
            RLV = new PS::F64vec **[S];
            qcos_glb = new PS::F64 **[S];
            qsin_glb = new PS::F64 **[S];
            for (int i = 0; i < S; ++i) {
                RLV[i] = new PS::F64vec *[S];
                qcos_glb[i] = new PS::F64 *[S];
                qsin_glb[i] = new PS::F64 *[S];
                for (int j = 0; j < S; ++j) {
                    RLV[i][j] = RLV_buf + (i * S + j) * S;
                    qcos_glb[i][j] = qcos_glb_buf + (i * S + j) * S;
                    qsin_glb[i][j] = qsin_glb_buf + (i * S + j) * S;
                }
            }
#pragma omp parallel for
            for (auto ix = -NIMG; ix <= NIMG; ix++) {
                for (auto iy = -NIMG; iy <= NIMG; iy++) {
                    for (auto iz = -NIMG; iz <= NIMG; iz++) {
                        RLV[ix + NIMG][iy + NIMG][iz + NIMG] = 0.0;
                        qcos_loc[ix + NIMG][iy + NIMG][iz + NIMG] = 0.0;
                        qsin_loc[ix + NIMG][iy + NIMG][iz + NIMG] = 0.0;
                        qcos_glb[ix + NIMG][iy + NIMG][iz + NIMG] = 0.0;
                        qsin_glb[ix + NIMG][iy + NIMG][iz + NIMG] = 0.0;
                        auto rlv_tmp = PS::F64vec(ix * gx, iy * gy, iz * gz);
                        RLV[ix + NIMG][iy + NIMG][iz + NIMG] = rlv_tmp;
                        for (PS::S32 i = 0; i < n; i++) {
                            const auto qi = ptcl[i].getCharge();
                            const auto posi = ptcl[i].pos;
                            qcos_loc[ix + NIMG][iy + NIMG][iz + NIMG] += qi * std::cos(rlv_tmp * posi);
                            qsin_loc[ix + NIMG][iy + NIMG][iz + NIMG] += qi * std::sin(rlv_tmp * posi);
                        }
                    }
                }
            }
            PS::Comm::getSum((double *)(**qcos_glb), (double *)(**qcos_loc), (2 * NIMG + 1) * (2 * NIMG + 1) * (2 * NIMG + 1));
            PS::Comm::getSum((double *)(**qsin_glb), (double *)(**qsin_loc), (2 * NIMG + 1) * (2 * NIMG + 1) * (2 * NIMG + 1));
            first_call = false;
        }

        static const PS::F64 inv_4alpha2 = 1.0 / (4.0 * alpha * alpha);
#pragma omp parallel for
        for (PS::S32 i = 0; i < n; i++) {
            const PS::F64vec posi = ptcl[i].pos;
            PS::F64vec grad = 0.0;
            PS::F64 pot = 0.0;
            for (int nx = -NIMG; nx <= NIMG; nx++) {
                for (int ny = -NIMG; ny <= NIMG; ny++) {
                    for (int nz = -NIMG; nz <= NIMG; nz++) {
                        if (nx == 0 && ny == 0 && nz == 0) continue;
                        const auto rlv_tmp = RLV[nx + NIMG][ny + NIMG][nz + NIMG];
                        const auto rlv_tmp2 = rlv_tmp * rlv_tmp;
                        const auto qcos_tmp = qcos_glb[nx + NIMG][ny + NIMG][nz + NIMG];
                        const auto qsin_tmp = qsin_glb[nx + NIMG][ny + NIMG][nz + NIMG];
                        const auto H = exp(-rlv_tmp2 * inv_4alpha2) / rlv_tmp2;
                        grad -= H * (sin(rlv_tmp * posi) * qcos_tmp - cos(rlv_tmp * posi) * qsin_tmp) * rlv_tmp;
                        pot += H * (cos(rlv_tmp * posi) * qcos_tmp + sin(rlv_tmp * posi) * qsin_tmp);
                    }
                }
            }
            ptcl[i].grad += 4.0 * pi / (len_unit_cell.x * len_unit_cell.y * len_unit_cell.z) * grad;
            ptcl[i].pot += 2.0 * pi / (len_unit_cell.x * len_unit_cell.y * len_unit_cell.z) * pot * 2.0;
            //  the last factor 2 is for 1/2 in the potential energy formula 1/2 * (q0*pot0 + q1*pot1 + ...) in the Ewald summation
        }
    }
};

template <typename Tpi, typename Tpj, typename Tforce>
static void calc_p0(Tpi *ep_i, const PS::S32 n_ip, Tpj *ep_j, const PS::S32 n_jp, Tforce *force) {
    for (PS::S32 i = 0; i < n_ip; i++) {
        PS::F64vec xi = ep_i[i].getPos();
        PS::F64vec grad = 0.0;  // nabla(pot)
        PS::F64 poti = 0.0;
        for (PS::S32 j = 0; j < n_jp; j++) {
            PS::F64vec rij = xi - ep_j[j].getPos();
            if (rij.x == 0.0 && rij.y == 0.0 && rij.z == 0.0) continue;
            PS::F64 mo = ep_j[j].getCharge();
            PS::F64 r2 = rij * rij;
            PS::F64 r_inv = 1.0 / sqrt(r2);
            PS::F64 r_inv2 = r_inv * r_inv;
            PS::F64 r_inv3 = r_inv2 * r_inv;
            poti += mo * r_inv;
            grad -= mo * r_inv3 * rij;
        }
        force[i].grad += grad;
        force[i].pot += poti;
    }
}

template <typename Tpi, typename Tpj, typename Tforce>
static void calc_p1(Tpi *ep_i, const PS::S32 n_ip, Tpj *ep_j, const PS::S32 n_jp, Tforce *force) {
    for (PS::S32 i = 0; i < n_ip; i++) {
        PS::F64vec xi = ep_i[i].getPos();
        PS::F64vec grad = 0.0;  // nabla(pot)
        PS::F64 poti = 0.0;
        for (PS::S32 j = 0; j < n_jp; j++) {
            PS::F64vec rij = xi - ep_j[j].getPos();
            if (rij.x == 0.0 && rij.y == 0.0 && rij.z == 0.0) continue;
            PS::F64 mo = ep_j[j].getMonopole();
            PS::F64vec di = ep_j[j].getDipole();
            PS::F64 r2 = rij * rij;
            PS::F64 r_inv = 1.0 / sqrt(r2);
            PS::F64 r_inv2 = r_inv * r_inv;
            PS::F64 r_inv3 = r_inv2 * r_inv;
            poti += mo * r_inv + di * rij * r_inv3;
            grad -= ((mo * r_inv3) * rij) + r_inv3 * ((3.0 * di * rij * r_inv2) * rij - di);
        }
        force[i].grad += grad;
        force[i].pot += poti;
    }
}

template <typename Tpi, typename Tpj, typename Tforce>
static void calc_p2(Tpi *ep_i, const PS::S32 n_ip, Tpj *ep_j, const PS::S32 n_jp, Tforce *force) {
    for (PS::S32 i = 0; i < n_ip; i++) {
        PS::F64vec xi = ep_i[i].getPos();
        PS::F64vec grad = 0.0;  // nabla(pot)
        PS::F64 poti = 0.0;
        for (PS::S32 j = 0; j < n_jp; j++) {
            PS::F64vec rij = xi - ep_j[j].getPos();
            if (rij.x == 0.0 && rij.y == 0.0 && rij.z == 0.0) continue;
            PS::F64 mo = ep_j[j].getMonopole();
            PS::F64vec di = ep_j[j].getDipole();
            PS::F64mat qd = ep_j[j].getQuadrupole();
            PS::F64 tr = qd.getTrace();
            PS::F64vec qr(qd.xx * rij.x + qd.xy * rij.y + qd.xz * rij.z, qd.xy * rij.x + qd.yy * rij.y + qd.yz * rij.z,
                          qd.xz * rij.x + qd.yz * rij.y + qd.zz * rij.z);
            PS::F64 rqr = qr * rij;
            PS::F64 r2 = rij * rij;
            PS::F64 r_inv = 1.0 / sqrt(r2);
            PS::F64 r_inv2 = r_inv * r_inv;
            PS::F64 r_inv3 = r_inv2 * r_inv;
            PS::F64 r_inv5 = r_inv3 * r_inv2;
            PS::F64 r_inv7 = r_inv5 * r_inv2;
            poti += mo * r_inv + di * rij * r_inv3 - (0.5 * tr - 1.5 * rqr * r_inv2) * r_inv3;
            grad -= ((mo * r_inv3) * rij) + r_inv3 * ((3.0 * di * rij * r_inv2) * rij - di) - 1.5 * tr * r_inv5 * rij - 3.0 * r_inv5 * qr +
                    7.5 * rqr * r_inv7 * rij;
        }
        force[i].grad += grad;
        force[i].pot += poti;
    }
}

// For PM
#if defined(USE_PM)
inline PS::F64 cube(const PS::F64 x) { return x * x * x; }
inline PS::F64 sqr(const PS::F64 x) { return x * x; }
inline PS::F64 gfactor_S2(const PS::F64 rad, const PS::F64 eps_pm) {
    PS::F64 R;
    PS::F64 g;
    PS::F64 S;
    R = 2.0 * rad / eps_pm;
    R = (R > 2.0) ? 2.0 : R;
    S = R - 1.0;
    S = (S > 0.0) ? S : 0.0;
    g = 1.0 + cube(R) * (-1.6 + sqr(R) * (1.6 + R * (-0.5 + R * (0.15 * R - 12.0 / 35.0)))) -
        cube(S) * cube(S) * (3.0 / 35.0 + R * (18.0 / 35.0 + 0.2 * R));
    return g;
}

template <typename Tpi, typename Tpj, typename Tforce>
class calc_pp_force_for_PM {
   public:
    void operator()(Tpi *ep_i, const PS::S32 ni, Tpj *ep_j, const PS::S32 nj, Tforce *ppforce) {
        for (PS::S32 i = 0; i < ni; i++) {
            for (PS::S32 j = 0; j < nj; j++) {
                PS::F64vec dr = ep_i[i].pos - ep_j[j].pos;
                if (dr.x == 0.0 && dr.y == 0.0 && dr.z == 0.0) continue;
                PS::F64 rsq = dr * dr;
                PS::F64 rad = sqrt(rsq);
                PS::F64 gfact = gfactor_S2(rad, 3.0 / SIZE_OF_MESH);
                PS::F64 rinv = 1.0 / rad;
                PS::F64 mrinv3 = ep_j[j].getCharge() * cube(rinv);
                ppforce[i].grad -= dr * gfact * mrinv3;
            }
        }
    }
};
#endif
