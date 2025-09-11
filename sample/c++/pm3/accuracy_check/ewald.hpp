#include <cmath>
#include <iostream>
#include <numbers>
#include <particle_simulator.hpp>

// pot  is the specific potential energy per unit charge. In the case of gravity, sign is reversed.
// grad is -nabla(pot)/qi (i.e. specific force) . In the case of gravity, the sign of pot is reversed, so specforce is -grad.

class Ewald {
   public:
    Ewald(const PS::F64 a, const PS::F64ort pu) : alpha_(a), pos_unit_cell_(pu) {
        len_unit_cell_ = pos_unit_cell_.high_ - pos_unit_cell_.low_;
        inv_len_unit_cell_ = PS::F64vec(1.0 / len_unit_cell_.x, 1.0 / len_unit_cell_.y, 1.0 / len_unit_cell_.z);
        volume_unit_cell_ = len_unit_cell_.x * len_unit_cell_.y * len_unit_cell_.z;
    }
    /*
    template <typename Tptcl>
    void setQtot(Tptcl *ptcl, const PS::S32 n) {
        q_tot_ = 0.0;
        #pragma omp parallel for reduction(+ : q_tot_)
        for (PS::S32 i = 0; i < n; i++) {
            q_tot_ += ptcl[i].getCharge();
        }
        //q_tot_ = PS::Comm::getSum(q_tot_loc);
        std::cout<<"q_tot_= "<<q_tot_<<std::endl;
    }
    */
    template <int NIMG, typename Tptcl>
    void calcForceSR(Tptcl *ptcl, const PS::S32 n) {
#pragma omp parallel for
        for (PS::S32 i = 0; i < n; i++) {
            PS::F64vec grad = 0.0;
            PS::F64 pot = 0.0;
            const PS::F64vec posi = ptcl[i].pos;
            for (int nx = -NIMG; nx <= NIMG; nx++) {
                for (int ny = -NIMG; ny <= NIMG; ny++) {
                    for (int nz = -NIMG; nz <= NIMG; nz++) {
                        PS::F64vec shift = PS::F64vec(nx * len_unit_cell_.x, ny * len_unit_cell_.y, nz * len_unit_cell_.z);
                        for (int j = 0; j < n; j++) {
                            if (nx == 0 && ny == 0 && nz == 0 && i == j) continue;
                            const PS::F64vec posj = ptcl[j].pos + shift;
                            const PS::F64vec posij = getMinimumDistance(posi - posj) + shift;
                            calcForceSRImpl(posij, ptcl[j].getCharge(), grad, pot);
                        }
                    }
                }
            }
            // ptcl[i].acc -= ptcl[i].getCharge() / ptcl[i].mass * grad;
#if defined(GRAVITY)
            ptcl[i].force -= grad;
#else
            // ptcl[i].force += grad;
            ptcl[i].force += (ptcl[i].getCharge() / ptcl[i].mass) * grad;
#endif
            // self energy // the last factor 2 is for 1/2 in the potential energy formula 1/2 * (q0*pot0 + q1*pot1 + ...) in the Ewald summation
            pot -= ptcl[i].getCharge() * alpha_ / std::sqrt(pi_) * 2.0;
            // charged system correction
            pot -= pi_ * q_tot_ / (alpha_ * alpha_ * volume_unit_cell_);
            ptcl[i].pot += pot;
        }
    }

    template <int NIMG, typename Tptcl>
    void calcForceLR(Tptcl ptcl[], const int n) {
        static const PS::F64 gx = 2.0 * pi_ / len_unit_cell_.x;
        static const PS::F64 gy = 2.0 * pi_ / len_unit_cell_.y;
        static const PS::F64 gz = 2.0 * pi_ / len_unit_cell_.z;
        static PS::F64vec RLV[2 * NIMG + 1][2 * NIMG + 1][2 * NIMG + 1];  // Reciprocal Lattice Vector
        static PS::F64 qcos[2 * NIMG + 1][2 * NIMG + 1][2 * NIMG + 1];    // qcos = sum (q * cos(RLV * pos))
        static PS::F64 qsin[2 * NIMG + 1][2 * NIMG + 1][2 * NIMG + 1];    // qsin = sum (q * sin(RLV * pos))
#pragma omp parallel for
        for (auto ix = -NIMG; ix <= NIMG; ix++) {
            for (auto iy = -NIMG; iy <= NIMG; iy++) {
                for (auto iz = -NIMG; iz <= NIMG; iz++) {
                    auto rlv_tmp = PS::F64vec(ix * gx, iy * gy, iz * gz);
                    RLV[ix + NIMG][iy + NIMG][iz + NIMG] = rlv_tmp;
                    for (PS::S32 i = 0; i < n; i++) {
                        const auto qi = ptcl[i].getCharge();
                        const auto posi = ptcl[i].pos;
                        qcos[ix + NIMG][iy + NIMG][iz + NIMG] += qi * std::cos(rlv_tmp * posi);
                        qsin[ix + NIMG][iy + NIMG][iz + NIMG] += qi * std::sin(rlv_tmp * posi);
                    }
                }
            }
        }
        const auto inv_4alpha2 = 1.0 / (4.0 * alpha_ * alpha_);
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
                        const auto qcos_tmp = qcos[nx + NIMG][ny + NIMG][nz + NIMG];
                        const auto qsin_tmp = qsin[nx + NIMG][ny + NIMG][nz + NIMG];
                        const auto H = exp(-rlv_tmp2 * inv_4alpha2) / rlv_tmp2;
                        grad += H * (sin(rlv_tmp * posi) * qcos_tmp - cos(rlv_tmp * posi) * qsin_tmp) * rlv_tmp;
                        pot += H * (cos(rlv_tmp * posi) * qcos_tmp + sin(rlv_tmp * posi) * qsin_tmp);
                    }
                }
            }
            // ptcl[i].acc -= 4.0 * pi_ / (len_unit_cell_.x * len_unit_cell_.y * len_unit_cell_.z) * ptcl[i].getCharge() / ptcl[i].mass * grad;
#if defined(GRAVITY)
            ptcl[i].force -= 4.0 * pi_ / (len_unit_cell_.x * len_unit_cell_.y * len_unit_cell_.z) * grad;
#else
            // ptcl[i].force += 4.0 * pi_ / (len_unit_cell_.x * len_unit_cell_.y * len_unit_cell_.z) * grad;
            ptcl[i].force += 4.0 * pi_ / (len_unit_cell_.x * len_unit_cell_.y * len_unit_cell_.z) * (ptcl[i].getCharge() / ptcl[i].mass) * grad;
#endif
            ptcl[i].pot += 2.0 * pi_ / (len_unit_cell_.x * len_unit_cell_.y * len_unit_cell_.z) * pot * 2.0;
            // the last factor 2 is for 1/2 in the potential energy formula 1/2 * (q0*pot0 + q1*pot1 + ...) in the Ewald summation
        }
    }

   private:
    PS::F64 alpha_;
    PS::F64ort pos_unit_cell_;
    PS::F64vec len_unit_cell_;
    PS::F64vec inv_len_unit_cell_;
    PS::F64 volume_unit_cell_;
    PS::F64 q_tot_ = 0.0;
    static constexpr PS::F64 pi_ = std::numbers::pi;
    static constexpr PS::F64 inv_sqrtpi_ = 1.0 / std::sqrt(pi_);

    // posij = posi - posj
    // fabs(posij) < len_unit_cell -> posij
    // otherwise, posij < 0 -> posij + len_unit_cell
    //            posij > 0 -> posij - len_unit_cell
    PS::F64vec getMinimumDistance(const PS::F64vec &posij) {
        PS::F64vec ret(0.0);
        // ret.x = std::abs(posij.x) < len_unit_cell_ ? posij.x : posij.x - copysign(len_unit_cell_.x, posij.x);
        // ret.y = std::abs(posij.x) < len_unit_cell_ ? posij.x : posij.x - copysign(len_unit_cell_.x, posij.x);
        // ret.z = std::abs(posij.x) < len_unit_cell_ ? posij.x : posij.x - copysign(len_unit_cell_.x, posij.x);
        ret.x = posij.x - std::round(posij.x * inv_len_unit_cell_.x) * len_unit_cell_.x;
        ret.y = posij.y - std::round(posij.y * inv_len_unit_cell_.y) * len_unit_cell_.y;
        ret.z = posij.z - std::round(posij.z * inv_len_unit_cell_.z) * len_unit_cell_.z;
        return ret;
    }
    // posij = posi - posj + Ln
    void calcForceSRImpl(const PS::F64vec &posij, const PS::F64 qj, PS::F64vec &grad, PS::F64 &pot) {
        auto dr2 = posij * posij;
        auto inv_dr = 1.0 / sqrt(dr2);
        auto inv_dr2 = inv_dr * inv_dr;
        auto dr = dr2 * inv_dr;
        auto alphadr = alpha_ * dr;
        grad += qj * (erfc(alphadr) * inv_dr + 2.0 * alpha_ * inv_sqrtpi_ * exp(-alphadr * alphadr)) * inv_dr2 * posij;
        pot += qj * erfc(alphadr) * inv_dr;
    }
};
