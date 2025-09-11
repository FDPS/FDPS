#pragma once

#include <iostream>
#include <iomanip>
#include "vector3.hpp"
#include "ps_defs.hpp"

namespace ParticleSimulator {

/*
  MSB is always 0.
  next 2bits represent quad index of 4 cells with level 1
 */
template <typename T>
class MortonKey {};

template <>
class MortonKey<TagMortonKeyNormal> {
   private:
    enum {
        kLevMax = 30,
    };
    F64 half_len_;
    F64vec center_;
    F64 normalized_factor_;
    U64 separateBit(const U64 _s_in) const {
        U64 _s = _s_in;
        _s = (_s | _s << 32) & 0x00000000ffffffff;   // 0000 0000 0000 0000 0000 0000 0000 0000 1111 1111 1111 1111 1111 1111 1111 1111
        _s = (_s | _s << 16) & 0x0000ffff0000ffff;   // 0000 0000 0000 0000 1111 1111 1111 1111 0000 0000 0000 0000 1111 1111 1111 1111
        _s = (_s | _s << 8) & 0x00ff00ff00ff00ff;    // 0000 0000 1111 1111 0000 0000 1111 1111
        _s = (_s | _s << 4) & 0x0f0f0f0f0f0f0f0f;    // 0000 1111 0000 1111 0000 1111
        _s = (_s | _s << 2) & 0x3333333333333333;    // 00 11 00 11 00 11
        return (_s | _s << 1) & 0x5555555555555555;  // 0101 0101
    }

   public:
    // MortonKey(){};
    //~MortonKey(){};
    // MortonKey(const MortonKey &){};
    // MortonKey & operator = (const MortonKey &);
    void initialize(const F64 half_len, const F64vec& center = 0.0) {
        half_len_ = half_len;
        center_ = center;
        normalized_factor_ = (1.0 / (half_len * 2.0)) * (1 << kLevMax);
    }

    KeyT getKey(const F64vec& pos) const {
        const F64vec cen = center_;
        const F64 hlen = half_len_;
        const F64 nfactor = normalized_factor_;
        const U64 nx = (U64)((pos.x - cen.x + hlen) * nfactor);
        const U64 ny = (U64)((pos.y - cen.y + hlen) * nfactor);
        return (separateBit(nx) << 1 | separateBit(ny));
    }
    template <typename Tkey>
    S32 getCellID(const S32 lev, const Tkey mkey) const {
        const auto s = mkey >> ((kLevMax - lev) * 2);
        return (S32)((s.hi_ & 0x3));
    }
    F64ort getCorrespondingTreeCell(const F64ort& box) const {
        const auto low = getKey(box.low_);
        const auto high = getKey(box.high_);
        const auto tmp = high ^ low;
        S32 lev = 0;
        for (S32 i = TREE_LEVEL_LIMIT - 1; i >= 0; i--) {
            if ((((tmp >> i * 2).hi_) & 0x3) != 0) {
                lev = TREE_LEVEL_LIMIT - i - 1;
                break;
            }
        }
        F64ort ret;
        const F64vec cen = center_;
        const F64 hlen = half_len_;
        if (lev == 0) {
            ret.low_.x = cen.x - hlen;
            ret.low_.y = cen.y - hlen;
            ret.high_.x = cen.x + hlen;
            ret.high_.y = cen.y + hlen;
        } else {
            F64vec cen_new = cen;
            F64 hlen_new = hlen;
            for (S32 i = 0; i < lev; i++) {
                S32 id = getCellID(i + 1, low);
                cen_new += hlen_new * SHIFT_CENTER[id];
                hlen_new *= 0.5;
            }
            ret.low_.x = cen_new.x - hlen_new;
            ret.low_.y = cen_new.y - hlen_new;
            ret.high_.x = cen_new.x + hlen_new;
            ret.high_.y = cen_new.y + hlen_new;
        }
        return ret;
    }
    F64 getCorrespondingFullLength(const S32 lev) const {
        F64 len = half_len_ * 2.0;
        for (S32 i = 0; i < lev; i++) {
            len *= 0.5;
        }
        return len;
    }
    template <typename Tkey>
    F64ort getPosTreeCell(const S32 lev, const Tkey& key) const {
        F64vec cen = center_;
        F64 hlen = half_len_;
        for (S32 i = 0; i < lev; i++) {
            S32 id = getCellID(i + 1, key);
            cen += hlen * SHIFT_CENTER[id];
            hlen *= 0.5;
        }
        F64ort ret;
        ret.low_.x = cen.x - hlen;
        ret.low_.y = cen.y - hlen;
        ret.high_.x = cen.x + hlen;
        ret.high_.y = cen.y + hlen;
        return ret;
    }
};

template <>
class MortonKey<TagMortonKeyMeshBased> {
   private:
    enum {
        kLevMax = 30,
    };
    F64vec half_len_;
    F64vec center_;
    F64vec pos_ref_;
    F64vec wid_ref_;
    F64vec wid_btm_;
    S32vec idx_ref_;
    S32 lev_ref_;
    U64 separateBit(const U64 _s_in) const {
        U64 _s = _s_in;
        _s = (_s | _s << 32) & 0x00000000ffffffff;   // 0000 0000 0000 0000 0000 0000 0000 0000 1111 1111 1111 1111 1111 1111 1111 1111
        _s = (_s | _s << 16) & 0x0000ffff0000ffff;   // 0000 0000 0000 0000 1111 1111 1111 1111 0000 0000 0000 0000 1111 1111 1111 1111
        _s = (_s | _s << 8) & 0x00ff00ff00ff00ff;    // 0000 0000 1111 1111 0000 0000 1111 1111
        _s = (_s | _s << 4) & 0x0f0f0f0f0f0f0f0f;    // 0000 1111 0000 1111 0000 1111
        _s = (_s | _s << 2) & 0x3333333333333333;    // 00 11 00 11 00 11
        return (_s | _s << 1) & 0x5555555555555555;  // 0101 0101
    }
    U64 shrinkBit(const U64 _s_in) const {
        U64 _s = _s_in;
        _s = (_s | (_s >> 1)) & 0x3333333333333333ULL;
        _s = (_s | (_s >> 2)) & 0x0f0f0f0f0f0f0f0fULL;
        _s = (_s | (_s >> 4)) & 0x00ff00ff00ff00ffULL;
        _s = (_s | (_s >> 8)) & 0x0000ffff0000ffffULL;
        _s = (_s | (_s >> 16)) & 0x00000000ffffffffULL;
        return _s;
    }

   public:
    void initialize(const F64 half_len, const F64vec& center) { /* Do nothing */ }
    void initialize(const F64vec half_len, const F64vec& center, const F64vec& pos_ref, const S32vec& idx_ref, const S32 lev_ref,
                    const F64vec& wid_ref) {
        half_len_ = half_len;
        center_ = center;
        pos_ref_ = pos_ref;
        idx_ref_ = idx_ref;
        wid_ref_ = wid_ref;
        lev_ref_ = lev_ref;
        wid_btm_ = wid_ref / (double)(1UL << (kLevMax - lev_ref));
    }

    KeyT getKey(const F64vec& pos) const {
        // Compute cell index
        S32vec idx_pm = GetCellIDMeasuredInUnitCell(pos_ref_, wid_ref_, pos);
        S32vec idx_btm;
        if (kLevMax > lev_ref) {
            // Compute sub-cell index in cell idx_pm.
            F64vec pos_cell;
            pos_cell.x = pos_ref_.x + wid_ref_.x * idx_pm.x;
            pos_cell.y = pos_ref_.y + wid_ref_.y * idx_pm.y;
            // [Notes]
            //    In principal, we can choose pos_cell and wid_btm
            //    different from those we are currently adopting
            //    because it is guaranteed that the particle of interest
            //    belongs to the cell to which it should belong.
            //    The only condition imposed on sub-cell indices is that
            //    the order of the particles belonging to that cell is
            //    preserved.
            idx_btm = GetCellIDMeasuredInUnitCell(pos_cell, wid_btm_, pos);
            const S32 maxval = (1 << (kLevMax - lev_ref_)) - 1;
            if (idx_btm.x < 0) idx_btm.x = 0;
            if (idx_btm.x > maxval) idx_btm.x = maxval;
            if (idx_btm.y < 0) idx_btm.y = 0;
            if (idx_btm.y > maxval) idx_btm.y = maxval;
        } else {
            idx_btm.x = idx_btm.y = 0;
        }
        U64 nx = (U64)((idx_pm.x + idx_ref_.x) * (1 << (kLevMax - lev_ref_)) + idx_btm.x);
        U64 ny = (U64)((idx_pm.y + idx_ref_.y) * (1 << (kLevMax - lev_ref_)) + idx_btm.y);

        // avoid to overflow
        nx = (((nx >> 63) & 0x1) == 1) ? 0 : nx;
        ny = (((ny >> 63) & 0x1) == 1) ? 0 : ny;
        U64 nmax = (1UL << kLevMax) - 1;
        nx = (nx > nmax) ? nmax : nx;
        ny = (ny > nmax) ? nmax : ny;

        KeyT ret;
        ret.hi_ = (separateBit(nx) << 1 | separateBit(ny));
        return ret;
    }
    template <typename Tkey>
    S32 getCellID(const S32 lev, const Tkey mkey) const {
        const auto s = mkey >> ((kLevMax - lev) * 2);
        return (S32)((s.hi_ & 0x3));
    }
    F64ort getCorrespondingTreeCell(const F64ort& box) const {
        const auto low = getKey(box.low_);
        const auto high = getKey(box.high_);
        const auto tmp = high ^ low;
        S32 lev = 0;
        for (S32 i = TREE_LEVEL_LIMIT - 1; i >= 0; i--) {
            if ((((tmp >> i * 2).hi_) & 0x3) != 0) {
                lev = TREE_LEVEL_LIMIT - i - 1;
                break;
            }
        }
        F64ort ret;
        const F64vec cen = center_;
        const F64vec hlen = half_len_;
        if (lev == 0) {
            ret.low_.x = cen.x - hlen.x;
            ret.low_.y = cen.y - hlen.y;
            ret.high_.x = cen.x + hlen.x;
            ret.high_.y = cen.y + hlen.y;
        } else {
            F64vec cen_new = cen;
            F64vec hlen_new = hlen;
            for (S32 i = 0; i < lev; i++) {
                S32 id = getCellID(i + 1, low);
                cen_new += hlen_new * SHIFT_CENTER[id];
                hlen_new *= 0.5;
            }
            ret.low_.x = cen_new.x - hlen_new.x;
            ret.low_.y = cen_new.y - hlen_new.y;
            ret.high_.x = cen_new.x + hlen_new.x;
            ret.high_.y = cen_new.y + hlen_new.y;
        }
        return ret;
    }
    F64 getCorrespondingFullLength(const S32 lev) const {
        F64vec len = half_len_ * 2.0;
        for (S32 i = 0; i < lev; i++) {
            len *= 0.5;
        }
        return std::max(len.x, len.y);
    }
    template <typename Tkey>
    F64ort getPosTreeCell(const S32 lev, const Tkey& key) const {
        F64vec cen = center_;
        F64vec hlen = half_len_;
        for (S32 i = 0; i < lev; i++) {
            S32 id = getCellID(i + 1, key);
            cen += hlen * SHIFT_CENTER[id];
            hlen *= 0.5;
        }
        F64ort ret;
        ret.low_.x = cen.x - hlen.x;
        ret.low_.y = cen.y - hlen.y;
        ret.high_.x = cen.x + hlen.x;
        ret.high_.y = cen.y + hlen.y;
        return ret;
    }
    // necessary functions for the PMMM feature
    S32vec getPMCellID(const F64vec pos) const {
        const KeyT mkey = getKey(pos);
        return getPMCellID(mkey);
    }
    template <typename Tkey>
    S32vec getPMCellID(const Tkey mkey) const {
        const U64 x = (mkey & 0xaaaaaaaaaaaaaaaaULL) >> 1;  // Repeat of pattern 10
        const U64 y = mkey & 0x5555555555555555ULL;         // Repeat of pattern 01
        S32 ix = (S32)(shrinkBit(x) >> (kLevMax - lev_ref_)) - idx_ref_.x;
        S32 iy = (S32)(shrinkBit(y) >> (kLevMax - lev_ref_)) - idx_ref_.y;
        return S32vec(ix, iy);
    }
};
}  // namespace ParticleSimulator
