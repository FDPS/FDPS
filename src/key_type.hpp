#pragma once

namespace ParticleSimulator {
#if defined(PARTICLE_SIMULATOR_USE_64BIT_KEY)
class KeyT {
   public:
    U64 hi_;
    //        KeyT() : hi_(0){}
    KeyT() {}
    KeyT(const U64 hi) : hi_(hi) {}
    void init() { hi_ = 0; }
    bool operator==(const KeyT &rhs) const { return (hi_ == rhs.hi_); }
    bool operator!=(const KeyT &rhs) const { return !(*this == rhs); }
    bool operator<(const KeyT &rhs) const { return (hi_ < rhs.hi_); }
    bool operator<=(const KeyT &rhs) const { return (hi_ <= rhs.hi_); }
    bool operator>(const KeyT &rhs) const { return !(*this <= rhs); }
    bool operator>=(const KeyT &rhs) const { return !(*this < rhs); }

    KeyT operator&(const KeyT &rhs) const { return KeyT(hi_ & rhs.hi_); }
    KeyT operator^(const KeyT &rhs) const { return KeyT(hi_ ^ rhs.hi_); }
    KeyT operator|(const KeyT &rhs) const { return KeyT(hi_ | rhs.hi_); }
    const KeyT &operator|=(const KeyT &rhs) {
        (this->hi_) |= rhs.hi_;
        return (*this);
    }
    // cast
    operator int() const { return (int)hi_; }
    operator unsigned int() const { return (unsigned int)hi_; }
    operator long() const { return (long)hi_; }
    operator unsigned long() const { return (unsigned long)hi_; }
    operator double() const { return (double)hi_; }

    template <typename T>
    KeyT operator<<(const T &s) const {
#if !defined(KEY_NEW_OPERATOR)
        assert(0);
#endif
        return KeyT(hi_ << s);
    }
    template <typename T>
    KeyT operator>>(const T &s) const {
#if !defined(KEY_NEW_OPERATOR)
        assert(0);
#endif
        return KeyT(hi_ >> s);
    }

    void dump(std::ostream &fout = std::cerr) const {
        fout << std::oct;
        fout << std::setw(TREE_LEVEL_LIMIT) << std::setfill('0') << std::oct << hi_ << std::endl;
        fout << std::dec;
    }

    friend std::ostream &operator<<(std::ostream &c, const KeyT &k) {
        // #if !defined(KEY_NEW_OPERATOR)
        //         assert(0);
        // #endif
        // c << k.hi_;
        c << std::oct;
        c << std::setw(TREE_LEVEL_LIMIT) << std::setfill('0') << k.hi_ << std::endl;
        c << std::dec;
        return c;
    }
};
#else  // for 96bit or 128bit
class KeyT {
   public:
    U64 hi_;
#ifdef PARTICLE_SIMULATOR_USE_96BIT_KEY
    U32 lo_;
    KeyT(const U64 hi, const U32 lo) : hi_(hi), lo_(lo) {}
#else
    U64 lo_;
    KeyT(const U64 hi, const U64 lo) : hi_(hi), lo_(lo) {}
#endif
    //        KeyT() : hi_(0), lo_(0){}
    KeyT() {}
    template <typename Tint>
    KeyT(const Tint lo) : hi_(0), lo_(lo) {}
    void init() {
        hi_ = 0;
        lo_ = 0;
    }
    bool operator==(const KeyT &rhs) const { return (hi_ == rhs.hi_) && (lo_ == rhs.lo_); }
    bool operator!=(const KeyT &rhs) const { return !(*this == rhs); }
    bool operator<(const KeyT &rhs) const { return (hi_ < rhs.hi_) || ((hi_ == rhs.hi_) && (lo_ < rhs.lo_)); }
    bool operator<=(const KeyT &rhs) const { return (hi_ < rhs.hi_) || ((hi_ == rhs.hi_) && (lo_ <= rhs.lo_)); }
    bool operator>(const KeyT &rhs) const { return !(*this <= rhs); }
    bool operator>=(const KeyT &rhs) const { return !(*this < rhs); }
    KeyT operator&(const KeyT &rhs) const { return KeyT(hi_ & rhs.hi_, lo_ & rhs.lo_); }
    KeyT operator^(const KeyT &rhs) const { return KeyT(hi_ ^ rhs.hi_, lo_ ^ rhs.lo_); }
    KeyT operator|(const KeyT &rhs) const { return KeyT(hi_ | rhs.hi_, lo_ | rhs.lo_); }
    const KeyT &operator|=(const KeyT &rhs) {
        (this->hi_) |= rhs.hi_;
        (this->lo_) |= rhs.lo_;
        return (*this);
    }
    // cast
    operator int() const { return (int)lo_; }
    operator unsigned int() const { return (unsigned int)lo_; }
    operator long() const { return (long)lo_; }
    operator unsigned long() const { return (unsigned long)lo_; }

    // shift
    template <typename T>
    KeyT operator<<(const T &s) const {
#if !defined(KEY_NEW_OPERATOR)
        assert(0);
#endif
#if defined(PARTICLE_SIMULATOR_USE_96BIT_KEY)
        U32 lo = ((U64)s < sizeof(lo_) * 8) ? lo_ << s : 0;
        U64 lo64 = (U64)lo_;
        const U64 rem = ((U64)s < sizeof(hi_) * 8) ? hi_ << s : 0;
        const U64 inc = ((U64)s < sizeof(lo_) * 8) ? (lo64 >> (sizeof(lo_) * 8 - s)) : lo64 << (s - sizeof(lo_) * 8);
        U64 hi = rem | inc;
#else
        U64 lo = ((U64)s < KEY_LEVEL_MAX_LO * 3) ? lo_ << s : 0;
        const U64 rem = (s == 0) ? 0 : (lo_ >> (KEY_LEVEL_MAX_LO * 3 - s));
        U64 hi = ((U64)s < KEY_LEVEL_MAX_LO * 3) ? ((hi_ << s) | rem) : lo_ << (s - KEY_LEVEL_MAX_LO * 3);
/*
            U64 lo = ((U64)s < sizeof(lo_)*8) ? lo_<<s : 0;
            const U64 rem = (s==0) ? 0 : (lo_>>(sizeof(hi_)*8-s));
            U64 hi = ((U64)s < sizeof(lo_)*8) ? ((hi_<<s) | rem) : lo_ << (s-sizeof(lo_)*8);
            */
#endif
        return KeyT(hi, lo);
    }
    template <typename T>
    KeyT operator>>(const T &s) const {
#if !defined(KEY_NEW_OPERATOR)
        assert(0);
#endif
#if defined(PARTICLE_SIMULATOR_USE_96BIT_KEY)
        U64 hi = ((U64)s < sizeof(hi_) * 8) ? hi_ >> s : 0;
        const U32 rem = ((U64)s < sizeof(lo_) * 8) ? lo_ >> s : 0;
        const U32 inc = (s == 0)                     ? 0
                        : ((U64)s < sizeof(hi_) * 8) ? (U32)((hi_ << (sizeof(hi_) * 8 - s)) >> (sizeof(lo_) * 8))
                                                     : (U32)(hi_ >> (s - sizeof(lo_) * 8));
        U32 lo = rem | inc;
#else
        // assert(s < 128);
        U64 hi = ((U64)s < KEY_LEVEL_MAX_HI * 3) ? hi_ >> s : 0;
        const U64 rem = (s == 0) ? 0 : (hi_ << (KEY_LEVEL_MAX_HI * 3 - s));
        U64 lo = ((U64)s < KEY_LEVEL_MAX_HI * 3) ? (lo_ >> s | rem) : hi_ >> (s - KEY_LEVEL_MAX_HI * 3);
        /*
        U64 hi = ((U64)s < sizeof(hi_)*8) ? hi_>>s : 0;
        const U64 rem = (s==0) ? 0 : (hi_<<(sizeof(hi_)*8-s));
        U64 lo = ((U64)s < sizeof(hi_)*8) ? (lo_>>s | rem) : hi_>>(s-sizeof(hi_)*8);
        */
#endif
        return KeyT(hi, lo);
    }

#if 1
    // get/set for external sort
    inline U64 get_hi_key() { return hi_; }
    inline U64 get_lo_key() { return lo_; }
    inline void set_hi_key(U64 val) { hi_ = val; }
    inline void set_lo_key(U64 val) { lo_ = val; }
#endif
    void dump(std::ostream &fout = std::cerr) const {
        U64 hi_tmp = hi_ & 0x7fffffffffffffff;
#if defined(PARTICLE_SIMULATOR_USE_96BIT_KEY)
        U32 lo_tmp = lo_ & 0x7fffffff;
#else
        U64 lo_tmp = lo_ & 0x7fffffffffffffff;
#endif
        fout << std::oct;
        fout << std::setw(KEY_LEVEL_MAX_HI) << std::setfill('0') << hi_tmp << " " << std::setw(KEY_LEVEL_MAX_LO) << std::setfill('0') << lo_tmp
             << std::endl;
        fout << std::dec;
    }

    friend std::ostream &operator<<(std::ostream &c, const KeyT &k) {
        /*
        c << std::oct;
        c << k.hi_ << "   " << k.lo_;
        c << std::dec;
        */
        U64 hi_tmp = k.hi_ & 0x7fffffffffffffff;
#if defined(PARTICLE_SIMULATOR_USE_96BIT_KEY)
        U32 lo_tmp = k.lo_ & 0x7fffffff;
#else
        U64 lo_tmp = k.lo_ & 0x7fffffffffffffff;
#endif
        c << std::oct;
        c << std::setw(KEY_LEVEL_MAX_HI) << std::setfill('0') << hi_tmp << " " << std::setw(KEY_LEVEL_MAX_LO) << std::setfill('0') << lo_tmp;
        c << std::dec;
        return c;
    }
};
#endif 
}  // namespace ParticleSimulator