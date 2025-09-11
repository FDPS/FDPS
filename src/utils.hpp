
namespace ParticleSimulator {

////////////////////////
// variadic templates //
/*
template <typename F, typename... Objs>
void InvokeOnEach(F&& f, Objs&&... objs) {
    (std::invoke(std::forward<F>(f), std::forward<Objs>(objs)), ...);
}


// Helper type to represent void-return results when collecting mixed return types
struct VoidResult {};

// invoke and normalize void / non-void results to be collectable in a tuple
template <class F, class Obj>
auto InvokeAndWrap(F&& f, Obj&& o) {
    using R = std::invoke_result_t<F, Obj>;
    if constexpr (std::is_void_v<R>) {
        std::invoke(std::forward<F>(f), std::forward<Obj>(o));
        return VoidResult{};
    } else {
        return std::invoke(std::forward<F>(f), std::forward<Obj>(o));
    }
}

// Collect results of invoking F on each object into a tuple. Void results become VoidResult.
template <class F, class... Objs>
auto CallAndCollect(F&& f, Objs&&... objs) {
    return std::tuple{InvokeAndWrap(std::forward<F>(f), std::forward<Objs>(objs))...};
}
*/



// ---------------------------------------------------------------------------

/*
// PS_INVOKE_ON_EACH(func, obj1, obj2, obj3);
#define PS_INVOKE_ON_EACH(METHOD, ...) \
    ParticleSimulator::InvokeOnEach([](auto &&o) { o.METHOD(); }, __VA_ARGS__)

//   PS_INVOKE_ON_EACH_WITH_METHOD_ARGS(func, args, obj1, obj2);
#define PS_INVOKE_ON_EACH_WITH_ARGS(METHOD, ARGS, ...) \
    ParticleSimulator::InvokeOnEach([&](auto &&o) { o.METHOD ARGS; }, __VA_ARGS__)

    
//   auto tup = PS_COLLECT_METHOD(getTimeProfile, obj1, obj2);    
#define PS_COLLECT_METHOD(METHOD, ...) \
    ParticleSimulator::CallAndCollect([](auto &o) { return o.METHOD(); }, __VA_ARGS__)

//   auto tup = PS_COLLECT_METHOD_ARGS(getValue, (arg1), objA, objB);
#define PS_COLLECT_METHOD_ARGS(METHOD, ARGS, ...) \
    ParticleSimulator::CallAndCollect([&](auto &o) { return o.METHOD ARGS; }, __VA_ARGS__)
*/


//////////////////////////
// small math functions //
template <typename Tfunc, typename T>
T MyReduction(Tfunc func, T &x0, T &x1) {
    return func(x0, x1);
}
template <typename Tfunc, typename T, typename... Args>
T MyReduction(Tfunc func, T &x, Args &...args) {
    return func(x, MyReduction(func, args...));
}
template <typename T, typename... Args>
T MyMax(T &x, Args &...args) {
    return MyReduction([](const auto &l, const auto &r) { return std::max(l, r); }, x, args...);
}
template <typename T, typename... Args>
T MyMin(T &x, Args &...args) {
    return MyReduction([](const auto &l, const auto &r) { return std::min(l, r); }, x, args...);
}
template <class T>
class Abs {
   public:
    T operator()(const T val) { return std::abs(val); }
};
static inline void Factorize(const S32 num, const S32 dim, S32 fact[]) {
    S32 num_tmp = num;
    std::vector<S32> num_vec(dim);
    for (S32 d = dim; d > 0; d--) {
        S32 tmp = (S32)pow((F64)num_tmp + 0.000001, (1.0 / d) * 1.000001);
        while (num_tmp % tmp != 0) {
            tmp--;
        }
        num_vec[d - 1] = tmp;
        num_tmp /= num_vec[d - 1];
    }
    std::sort(num_vec.begin(), num_vec.end(), std::greater<S32>());
    for (auto k = 0; k < dim; k++) {
        fact[k] = num_vec[k];
    }
    S32 num_check = fact[0];
    for (S32 d = 1; d < dim; d++) {
        num_check *= fact[d];
        assert(fact[d - 1] >= fact[d]);
    }
    assert(num_check == num);
}
template <typename T>
static inline void PrefixSum(const ReallocatableArray<T> &n, ReallocatableArray<T> &n_disp) {
    n_disp.resizeNoInitialize(n.size() + 1);
    n_disp[0] = 0;
    for (int i = 0; i < n.size(); i++) {
        n_disp[i + 1] = n_disp[i] + n[i];
    }
}
// return the minimum integer lev such that 2^lev >= in
static PS_INLINE S32 GetMinExpOfTwoGE(const S32 in) {  // GE := equal to or greater than
    assert(in > 0);
    S32 tmp = 1, lev = 0;
    while (tmp < in) {
        tmp *= 2;
        lev++;
    }
    return lev;
}

// return the minimum integer lev such that 2^lev > in
static PS_INLINE S32 GetMinExpOfTwoGT(const S32 in) {  // GT := greater than
    assert(in > 0);
    S32 tmp = 1, lev = 0;
    while (tmp <= in) {
        tmp *= 2;
        lev++;
    }
    return lev;
}

//////////////
// FANCTORs //
struct Cmpvec {
    F64 F64vec::*loc;
    Cmpvec(F64 F64vec::*_loc) : loc(_loc) {}
    bool operator()(const F64vec &lhs, const F64vec &rhs) { return (lhs.*loc < rhs.*loc); }
};

////////////////////
/// BOX HANDLING ///
inline bool IsInBox(const F64vec &pos, const F64vec &center, const F64 half_length, const F64 tolerance = 1e-6) {
    const F64 tol = -std::fabs(tolerance);
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
    return ((pos.x - (center.x - half_length)) >= tol) && (((center.x + half_length) - pos.x) >= tol) &&
           ((pos.y - (center.y - half_length)) >= tol) && (((center.y + half_length) - pos.y) >= tol);
#else
    return ((pos.x - (center.x - half_length)) >= tol) && (((center.x + half_length) - pos.x) >= tol) &&
           ((pos.y - (center.y - half_length)) >= tol) && (((center.y + half_length) - pos.y) >= tol) &&
           ((pos.z - (center.z - half_length)) >= tol) && (((center.z + half_length) - pos.z) >= tol);
#endif
}
template <class Tp>
inline F64ort GetMinBoxSingleThread(const Tp ptcl[], const S32 n) {
    F64ort pos_box(ptcl[0].getPos(), ptcl[0].getPos());
    for (S32 i = 1; i < n; i++) {
        pos_box.merge(ptcl[i].getPos());
    }
    return pos_box;
}

/////////////////////
/// DATA HANDLING ///
// paking data into a single array
// data_in[][] -> data_out[]: data_out[0] = data_in[0][0], data_out[1] = data_in[0][1], ...data_out[n] = data_in[1][0], ...
template <typename T>
inline void PackData(const ReallocatableArray<T> *data_in, const S32 n_buf, ReallocatableArray<T> &data_out, const S32 offset = 0) {
    // openmp version
    S32 size_data = offset;
    PS_OMP(omp parallel for reduction(+:size_data))
    for (S32 i = 0; i < n_buf; i++) {
        size_data += data_in[i].size();
    }
    data_out.resizeNoInitialize(size_data);
    PS_OMP(omp parallel for)
    for (S32 i = 0; i < n_buf; i++) {
        S32 head_adr = offset;
        for (S32 j = 0; j < i; j++) {
            head_adr += data_in[j].size();
        }
        for (S32 j = 0; j < data_in[i].size(); j++) {
            data_out[j + head_adr] = data_in[i][j];
        }
    }
}
// NOTE: This function should be used only in OpenMP scopes
template <typename T>
inline void PackDataInOmp(const ReallocatableArray<T> *data_in, ReallocatableArray<T> &data_out, const S32 nth, const S32 ith, const S32 offset = 0) {
    if (ith == 0) {
        S32 size_data = offset;
        for (S32 i = 0; i < nth; i++) {
            size_data += data_in[i].size();
        }
        data_out.resizeNoInitialize(size_data);
    }
    PS_OMP_BARRIER
    S32 head_adr = offset;
    for (S32 j = 0; j < ith; j++) {
        head_adr += data_in[j].size();
    }
    for (S32 j = 0; j < data_in[ith].size(); j++) {
        data_out[j + head_adr] = data_in[ith][j];
    }
}

// split data
// [head, end)
inline void CalcAdrToSplitData(S32 &head, S32 &end, const S32 i_th, const S32 n_div, const S32 n_tot) {
    head = (n_tot / n_div) * i_th + std::min(n_tot % n_div, i_th);
    end = (n_tot / n_div) * (i_th + 1) + std::min(n_tot % n_div, (i_th + 1));
}

inline std::tuple<S32, S32> CalcAdrToSplitData(const S32 i_th, const S32 n_div, const S32 n_tot) {
    S32 head, end;
    CalcAdrToSplitData(head, end, i_th, n_div, n_tot);
    return {head, end};
}




///////////////////
// MSB HANDLING ///
template <class T>
inline T GetMSB(const T val);
template <>
inline U64 GetMSB(const U64 val) {
    return (val >> 63) & 0x1;
}
template <>
inline U32 GetMSB(const U32 val) {
    return (val >> 31) & 0x1;
}

template <class T>
inline T ClearMSB(const T val);
template <>
inline U64 ClearMSB(const U64 val) {
    return val & 0x7fffffffffffffff;
}
template <>
inline U32 ClearMSB(const U32 val) {
    return val & 0x7fffffff;
}

template <class T>
inline T SetMSB(const T val);
template <>
inline U64 SetMSB(const U64 val) {
    return val | 0x8000000000000000;
}
template <>
inline U32 SetMSB(const U32 val) {
    return val | 0x80000000;
}

//////////////////////////////
/// UNIQUE RANDOM SAMPLING ///
/*
template <typename Tsrc, typename Tdst, typename Tcopy>
inline void RandomSamplingOMP(Tsrc *val_src, Tdst *val_dst, const S32 n_src, const S32 n_dst, Tcopy copy) {
    PS_OMP_PARALLEL {
        auto i_th = Comm::getThreadNum();
        auto n_div = Comm::getNumberOfThread();
        // Calculate the source address ranges
        auto [head_src, end_src] = CalcAdrToSplitData(i_th, n_div, n_src);
        // Calculate the destination address ranges
        auto [head_dst, end_dst] = CalcAdrToSplitData(i_th, n_div, n_dst);
        auto n_src_loc = end_src - head_src;
        auto n_dst_loc = end_dst - head_dst;
        RandomSampling<Tsrc, Tdst, Tcopy>(val_src + head_src, val_dst + head_dst, n_src_loc, n_dst_loc, copy);
    }
}
*/
template <typename Tsrc, typename Tdst, typename Tcopy>
inline void RandomSampling(Tsrc *val_src, Tdst *val_dst, const S32 n_src, const S32 n_dst, Tcopy copy) {
    assert(n_dst <= n_src);
    thread_local std::mt19937 mt(Comm::getRank() * Comm::getNumberOfThread() + Comm::getThreadNum());
    if (n_src == 0) return;
    std::uniform_int_distribution<S32> dist(0, n_src - 1);
    ReallocatableArray<S32> record(n_dst, n_dst, MemoryAllocMode::Pool);
    for (S32 i = 0; i < n_dst; i++) {
        S32 j = dist(mt);
        Tsrc hold = val_src[j];
        val_src[j] = val_src[i];
        val_src[i] = hold;
        record[i] = j;
    }
    for (S32 i = 0; i < n_dst; i++) {
        copy(val_src[i], val_dst[i]);
    }
    for (S32 i = n_dst - 1; i >= 0; i--) {
        S32 j = record[i];
        Tsrc hold = val_src[j];
        val_src[j] = val_src[i];
        val_src[i] = hold;
    }
}

template <typename Tsrc, typename Tdst, typename Tcopy>
inline void RandomSampling(Tsrc *val_src, ReallocatableArray<Tdst> &val_dst, const S32 n_src, Tcopy copy) {
    const auto n_dst = val_dst.size();
    assert(n_dst <= n_src);
    thread_local std::mt19937 mt(Comm::getRank() * Comm::getNumberOfThread() + Comm::getThreadNum());
    if (n_src == 0) return;
    std::uniform_int_distribution<S32> dist(0, n_src - 1);
    ReallocatableArray<S32> record(n_dst, n_dst, MemoryAllocMode::Pool);
    for (S32 i = 0; i < n_dst; i++) {
        S32 j = dist(mt);
        Tsrc hold = val_src[j];
        val_src[j] = val_src[i];
        val_src[i] = hold;
        record[i] = j;
    }
    for (S32 i = 0; i < n_dst; i++) {
        copy(val_src[i], val_dst[i]);
    }
    for (S32 i = n_dst - 1; i >= 0; i--) {
        S32 j = record[i];
        Tsrc hold = val_src[j];
        val_src[j] = val_src[i];
        val_src[i] = hold;
    }
}

template <typename Tsrc, typename Tdst, typename Tcopy>
inline void RandomSampling(ReallocatableArray<Tsrc> &val_src_ar, ReallocatableArray<Tdst> &val_dst_ar, Tcopy copy, const S32 offset_src = 0,
                           const S32 offset_dst = 0) {
    const auto n_src = val_src_ar.size() - offset_src;
    const auto n_dst = val_dst_ar.size() - offset_dst;
    assert(n_dst <= n_src);
    thread_local std::mt19937 mt(Comm::getRank() * Comm::getNumberOfThread() + Comm::getThreadNum());
    // thread_local std::mt19937 mt(0);
    if (n_src == 0) return;
    std::uniform_int_distribution<S32> dist(0, n_src - 1);
    ReallocatableArray<S32> record(n_dst, n_dst, MemoryAllocMode::Pool);
    Tsrc *val_src = val_src_ar.getPointer(offset_src);
    Tdst *val_dst = val_dst_ar.getPointer(offset_dst);
    for (S32 i = 0; i < n_dst; i++) {
        S32 j = dist(mt);
        Tsrc hold = val_src[j];
        val_src[j] = val_src[i];
        val_src[i] = hold;
        record[i] = j;
    }
    for (S32 i = 0; i < n_dst; i++) {
        copy(val_src[i], val_dst[i]);
    }
    for (S32 i = n_dst - 1; i >= 0; i--) {
        S32 j = record[i];
        Tsrc hold = val_src[j];
        val_src[j] = val_src[i];
        val_src[i] = hold;
    }
}



//////////////////////////////////////////////////////
/// subtraction for both periodic and non-periodic ///
template <typename T>
PS_INLINE T GetSubImpl(const T l, const T r, const T len_peri) {
    const T diff = l - r;
    const T d0 = std::abs(diff);
    const T d1 = len_peri - d0;
    const T ret = d0 < d1 ? diff : (diff < T(0) ? d1 : -d1);
    return ret;
}
template <typename T>
PS_INLINE T GetSubImpl(const T l, const T r) {
    const T diff = l - r;
    return diff;
}
template <int bits>
PS_INLINE F64 GetSub(const F64 l, const F64 r, const F64 len_peri) {
    Abort();
    return -1;
}
template <>
PS_INLINE F64 GetSub<0>(const F64 l, const F64 r, const F64 len_peri) {
    return GetSubImpl(l, r);
}
template <>
PS_INLINE F64 GetSub<1>(const F64 l, const F64 r, const F64 len_peri) {
    return GetSubImpl(l, r, len_peri);
}
template <int bits>
PS_INLINE F64vec GetSub(const F64vec l, const F64vec r, const F64vec len_peri) {
    auto dx = GetSub<((bits >> 0) & 0x1)>(l.x, r.x, len_peri.x);
    auto dy = GetSub<((bits >> 1) & 0x1)>(l.y, r.y, len_peri.y);
#if !defined PARTICLE_SIMULATOR_TWO_DIMENSION
    auto dz = GetSub<((bits >> 2) & 0x1)>(l.z, r.z, len_peri.z);
    return F64vec(dx, dy, dz);
#else
    return F64vec(dx, dy);
#endif
}
/// subtraction for both periodic and non-periodic ///
//////////////////////////////////////////////////////

///////////////////////
///   calc distance ///
// between points
template <typename T>
PS_INLINE T GetDistanceMin1DImpl(const T a, const T b) {
    return std::abs(a - b);
}
template <typename T>
PS_INLINE T GetDistanceMin1DPeriImpl(const T a, const T b, const T len) {
    const auto h0 = std::abs(a - b);
    const auto h1 = len - h0;
    return std::min(h0, h1);
}
template <int bit>
PS_INLINE F64 GetDistanceMin1D(const F64 a, const F64 b, const F64 len) {
    Abort();
    return -1;
}
template <>
PS_INLINE F64 GetDistanceMin1D<0>(const F64 a, const F64 b, const F64 len) {
    return GetDistanceMin1DImpl(a, b);
}

template <>
PS_INLINE F64 GetDistanceMin1D<1>(const F64 a, const F64 b, const F64 len) {
    return GetDistanceMin1DPeriImpl(a, b, len);
}

//////////////////
// between lines
PS_INLINE F64 GetDistanceMin1DImpl(const F64 al, const F64 ah, const F64 bl, const F64 bh) {
    const auto dx0 = al - bh;
    const auto dx1 = bl - ah;
    const auto dx = (dx0 * dx1 >= 0.0) ? 0.0 : ((dx0 > 0.0) ? dx0 : dx1);
    return dx;
}

PS_INLINE F64 GetDistanceMin1DPeriImpl(const F64 al, const F64 ah, const F64 bl, const F64 bh, const F64 len) {
    const auto dx0 = al - bh;
    const auto dx1 = bl - ah;
    const auto dis0 = std::max(dx0, dx1);
    const auto dis1 = std::max((len + std::min(dx0, dx1)), 0.0);
    const auto dis = (dx0 * dx1 >= 0.0) ? 0.0 : std::min(dis0, dis1);
    return dis;
}

template <int bit>
PS_INLINE F64 GetDistanceMin1D(const F64 al, const F64 ah, const F64 bl, const F64 bh, const F64 len) {
    Abort();
    return -1;
}
template <>
PS_INLINE F64 GetDistanceMin1D<0>(const F64 al, const F64 ah, const F64 bl, const F64 bh, const F64 len) {
    return GetDistanceMin1DImpl(al, ah, bl, bh);
}
template <>
PS_INLINE F64 GetDistanceMin1D<1>(const F64 al, const F64 ah, const F64 bl, const F64 bh, const F64 len) {
    return GetDistanceMin1DPeriImpl(al, ah, bl, bh, len);
}

//////////////////
// between point and line
PS_INLINE F64 GetDistanceMin1DImpl(const F64 al, const F64 ah, const F64 b) {
    const auto dx0 = b - ah;
    const auto dx1 = al - b;
    const auto dx = (dx0 * dx1 >= 0.0) ? 0.0 : ((dx0 > 0.0) ? dx0 : dx1);
    return dx;
}

PS_INLINE F64 GetDistanceMin1DPeriImpl(const F64 al, const F64 ah, const F64 b, const F64 len) {
    const auto dx0 = b - ah;
    const auto dx1 = al - b;
    const auto dis0 = std::max(dx0, dx1);
    // const auto dis1 = len + std::min(dx0, dx1);
    const auto dis1 = std::max((len + std::min(dx0, dx1)), 0.0);
    const auto dis = (dx0 * dx1 >= 0.0) ? 0.0 : std::min(dis0, dis1);
    return dis;
}

template <int bit>
PS_INLINE F64 GetDistanceMin1D(const F64 al, const F64 ah, const F64 b, const F64 len) {
    Abort();
    return -1;
}
template <>
PS_INLINE F64 GetDistanceMin1D<0>(const F64 al, const F64 ah, const F64 b, const F64 len) {
    return GetDistanceMin1DImpl(al, ah, b);
}
template <>
PS_INLINE F64 GetDistanceMin1D<1>(const F64 al, const F64 ah, const F64 b, const F64 len) {
    return GetDistanceMin1DPeriImpl(al, ah, b, len);
}

//////////////////
// between F64vec
PS_INLINE F64 GetDistanceMinSq(const F64vec &pos0, const F64vec &pos1) {
    const auto dpos = pos0 - pos1;
    return dpos * dpos;
}

template <int bits>
PS_INLINE F64 GetDistanceMinSq(const F64vec &pos0, const F64vec &pos1, const F64vec &len_peri) {
    const auto dx = GetDistanceMin1D<((bits >> 0) & 0x1)>(pos0.x, pos1.x, len_peri.x);
    const auto dy = GetDistanceMin1D<((bits >> 1) & 0x1)>(pos0.y, pos1.y, len_peri.y);
    auto dis = dx * dx + dy * dy;
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
    const auto dz = GetDistanceMin1D<((bits >> 2) & 0x1)>(pos0.z, pos1.z, len_peri.z);
    dis += dz * dz;
#endif
    return dis;
}

//////////////////
// between F64ort
PS_INLINE F64 GetDistanceMinSq(const F64ort &pos0, const F64ort &pos1) {
    const auto dx = GetDistanceMin1DImpl(pos0.low_.x, pos0.high_.x, pos1.low_.x, pos1.high_.x);
    const auto dy = GetDistanceMin1DImpl(pos0.low_.y, pos0.high_.y, pos1.low_.y, pos1.high_.y);
    auto dis = dx * dx + dy * dy;
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
    const auto dz = GetDistanceMin1DImpl(pos0.low_.z, pos0.high_.z, pos1.low_.z, pos1.high_.z);
    dis += dz * dz;
#endif
    return dis;
}

template <int bits = 0>
PS_INLINE F64 GetDistanceMinSq(const F64ort &pos0, const F64ort &pos1, const F64vec &len_peri) {
    const auto dx = GetDistanceMin1D<((bits >> 0) & 0x1)>(pos0.low_.x, pos0.high_.x, pos1.low_.x, pos1.high_.x, len_peri.x);
    const auto dy = GetDistanceMin1D<((bits >> 1) & 0x1)>(pos0.low_.y, pos0.high_.y, pos1.low_.y, pos1.high_.y, len_peri.y);
    auto dis = dx * dx + dy * dy;
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
    const auto dz = GetDistanceMin1D<((bits >> 2) & 0x1)>(pos0.low_.z, pos0.high_.z, pos1.low_.z, pos1.high_.z, len_peri.z);
    dis += dz * dz;
#endif
    return dis;
}

//////////////////
// between F64ort and F64vec
PS_INLINE F64 GetDistanceMinSq(const F64ort &pos0, const F64vec &pos1) {
    const auto dx = GetDistanceMin1DImpl(pos0.low_.x, pos0.high_.x, pos1.x);
    const auto dy = GetDistanceMin1DImpl(pos0.low_.y, pos0.high_.y, pos1.y);
    auto dis = dx * dx + dy * dy;
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
    const F64 dz = GetDistanceMin1DImpl(pos0.low_.z, pos0.high_.z, pos1.z);
    dis += dz * dz;
#endif
    return dis;
}

PS_INLINE F64 GetDistanceMinSq(const F64vec &pos0, const F64ort &pos1) { return GetDistanceMinSq(pos1, pos0); }

template <int bits>
PS_INLINE F64 GetDistanceMinSq(const F64ort &pos0, const F64vec &pos1, const F64vec &len_peri) {
    const auto dx = GetDistanceMin1D<((bits >> 0) & 0x1)>(pos0.low_.x, pos0.high_.x, pos1.x, len_peri.x);
    const auto dy = GetDistanceMin1D<((bits >> 1) & 0x1)>(pos0.low_.y, pos0.high_.y, pos1.y, len_peri.y);
    auto dis = dx * dx + dy * dy;
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
    const F64 dz = GetDistanceMin1D<((bits >> 2) & 0x1)>(pos0.low_.z, pos0.high_.z, pos1.z, len_peri.z);
    dis += dz * dz;
#endif
    return dis;
}

template <int bits>
PS_INLINE F64 GetDistanceMinSq(const F64vec &pos0, const F64ort &pos1, const F64vec &len_peri) {
    return GetDistanceMinSq<bits>(pos1, pos0, len_peri);
}

template <int bits>
PS_INLINE F64 GetPosNearestImage1D(const F64 &pos_target, const F64 &pos_org, const F64 &len_peri) {
    Abort();
    return -1;
}
template <>
PS_INLINE F64 GetPosNearestImage1D<0>(const F64 &pos_target, const F64 &pos_org, const F64 &len_peri) {
    return pos_target;
}
template <>
PS_INLINE F64 GetPosNearestImage1D<1>(const F64 &pos_target, const F64 &pos_org, const F64 &len_peri) {
    auto dx = pos_target - pos_org;
    auto dx_abs = std::abs(dx);
    auto pos_new = (dx_abs <= len_peri - dx_abs) ? pos_target : (dx > 0.0 ? pos_target - len_peri : pos_target + len_peri);
    return pos_new;
}
template <int bits>
PS_INLINE F64vec GetPosNearestImage(const F64vec &pos_target, const F64vec &pos_org, const F64vec &len_peri) {
    auto x_new = GetPosNearestImage1D<((bits >> 0) & 0x1)>(pos_target.x, pos_org.x, len_peri.x);
    auto y_new = GetPosNearestImage1D<((bits >> 1) & 0x1)>(pos_target.y, pos_org.y, len_peri.y);
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
    auto z_new = GetPosNearestImage1D<((bits >> 2) & 0x1)>(pos_target.z, pos_org.z, len_peri.z);
    return F64vec(x_new, y_new, z_new);
#else
    return F64vec(x_new, y_new);
#endif
}

//////////////////
// over lapped functions
#if __cplusplus >= PS_CPLUSPLUS_17
template <typename...>
constexpr bool FALSE_V = false;
template <typename T0, typename T1>
PS_INLINE bool IsOverlapped(const T0 &pos0, const T1 &pos1) {
    if constexpr ((std::is_same_v<T0, F64vec> && std::is_same_v<T1, F64vec>) || (std::is_same_v<T0, F64ort> && std::is_same_v<T1, F64ort>) ||
                  (std::is_same_v<T0, F64vec> && std::is_same_v<T1, F64ort>) || (std::is_same_v<T0, F64ort> && std::is_same_v<T1, F64vec>)) {
        // asm("#IF TRUE");
        return GetDistanceMinSq(pos0, pos1) <= 0.0;
    } else {
        // asm("#IF FALSE");
        static_assert(FALSE_V<T0, T1>);
    }
}
template <int bits, typename T0, typename T1>
PS_INLINE bool IsOverlapped(const T0 &pos0, const T1 &pos1, const F64vec &len_peri) {
    if constexpr ((std::is_same_v<T0, F64vec> && std::is_same_v<T1, F64vec>) || (std::is_same_v<T0, F64ort> && std::is_same_v<T1, F64ort>) ||
                  (std::is_same_v<T0, F64vec> && std::is_same_v<T1, F64ort>) || (std::is_same_v<T0, F64ort> && std::is_same_v<T1, F64vec>)) {
        return GetDistanceMinSq<bits>(pos0, pos1, len_peri) <= 0.0;
    } else {
        static_assert(FALSE_V<T0, T1>);
    }
}
#else
PS_INLINE bool IsOverlapped(const F64ort &pos0, const F64ort &pos1) { return GetDistanceMinSq(pos0, pos1) <= 0.0; }
PS_INLINE bool IsOverlapped(const F64vec &pos0, const F64vec &pos1) { return GetDistanceMinSq(pos0, pos1) <= 0.0; }
PS_INLINE bool IsOverlapped(const F64ort &pos0, const F64vec &pos1) { return GetDistanceMinSq(pos0, pos1) <= 0.0; }
PS_INLINE bool IsOverlapped(const F64vec &pos0, const F64ort &pos1) { return GetDistanceMinSq(pos0, pos1) <= 0.0; }

template <int bits>
PS_INLINE bool IsOverlapped(const F64ort &pos0, const F64ort &pos1, const F64vec &len_peri) {
    return GetDistanceMinSq<bits>(pos0, pos1, len_peri) <= 0.0;
}
template <int bits>
PS_INLINE bool IsOverlapped(const F64vec &pos0, const F64vec &pos1, const F64vec &len_peri) {
    return GetDistanceMinSq<bits>(pos0, pos1, len_peri) <= 0.0;
}
template <int bits>
PS_INLINE bool IsOverlapped(const F64ort &pos0, const F64vec &pos1, const F64vec &len_peri) {
    return GetDistanceMinSq<bits>(pos0, pos1, len_peri) <= 0.0;
}
template <int bits>
PS_INLINE bool IsOverlapped(const F64vec &pos0, const F64ort &pos1, const F64vec &len_peri) {
    return GetDistanceMinSq<bits>(pos0, pos1, len_peri) <= 0.0;
}
#endif

// for search modes other than SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE
static inline void CalcNumberAndShiftOfImageDomain(ReallocatableArray<F64vec> &shift_image_domain, const F64vec &size_root_domain,
                                                   const F64ort &pos_my_domain, const F64ort &pos_target_domain, const bool periodic_axis[],
                                                   const bool clear = true) {
    if (clear) shift_image_domain.clearSize();
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
    if (periodic_axis[0] == false && periodic_axis[1] == false) {
        shift_image_domain.push_back(F64vec(0.0));  // NOTE: sign is plus
        return;
    }
#else
    if (periodic_axis[0] == false && periodic_axis[1] == false && periodic_axis[2] == false) {
        shift_image_domain.push_back(F64vec(0.0));  // NOTE: sign is plus
        return;
    }
#endif
    if (pos_my_domain.overlapped(pos_target_domain)) {
        shift_image_domain.push_back(F64vec(0.0));  // NOTE: sign is plus
    }
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
    for (S32 lev = 1;; lev++) {
        S32 n_image_per_level = 0;
        S32 lev_x = 0;
        S32 lev_y = 0;
        if (periodic_axis[0] == true) lev_x = lev;
        if (periodic_axis[1] == true) lev_y = lev;
        for (S32 ix = -lev_x; ix <= lev_x; ix++) {
            for (S32 iy = -lev_y; iy <= lev_y; iy++) {
                if ((std::abs(ix) != lev_x || periodic_axis[0] == false) && (std::abs(iy) != lev_y || periodic_axis[1] == false)) continue;
                const F64vec shift_tmp(ix * size_root_domain.x, iy * size_root_domain.y);
                const F64ort pos_image_tmp = pos_target_domain.shift(shift_tmp);
                if (pos_my_domain.overlapped(pos_image_tmp)) {
                    shift_image_domain.push_back(shift_tmp);  // NOTE: sign is plus
                    n_image_per_level++;
                }
            }
        }
        if (n_image_per_level == 0) break;
    }
#else
    for (S32 lev = 1;; lev++) {
        S32 n_image_per_level = 0;
        S32 lev_x = 0;
        S32 lev_y = 0;
        S32 lev_z = 0;
        if (periodic_axis[0] == true) lev_x = lev;
        if (periodic_axis[1] == true) lev_y = lev;
        if (periodic_axis[2] == true) lev_z = lev;
        for (S32 ix = -lev_x; ix <= lev_x; ix++) {
            for (S32 iy = -lev_y; iy <= lev_y; iy++) {
                for (S32 iz = -lev_z; iz <= lev_z; iz++) {
                    if ((std::abs(ix) != lev_x || periodic_axis[0] == false) && (std::abs(iy) != lev_y || periodic_axis[1] == false) &&
                        (std::abs(iz) != lev_z || periodic_axis[2] == false))
                        continue;
                    const F64vec shift_tmp(ix * size_root_domain.x, iy * size_root_domain.y, iz * size_root_domain.z);
                    const F64ort pos_image_tmp = pos_target_domain.shift(shift_tmp);
                    if (pos_my_domain.overlapped(pos_image_tmp)) {
                        shift_image_domain.push_back(shift_tmp);  // NOTE: sign is plus
                        n_image_per_level++;
                    }
                }
            }
        }
        if (n_image_per_level == 0) break;
    }
#endif
}

// for SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE
template <typename Tmorton_key>
inline void CalcNumberAndShiftOfImageDomain(ReallocatableArray<F64vec> &shift_image_domain, const F64vec &size_root_domain,
                                            const F64ort &pos_my_domain, const F64ort &pos_target_domain, const bool periodic_axis[], const S32 icut,
                                            const Tmorton_key &morton_key, const bool clear = true) {
    if (clear) shift_image_domain.clearSize();
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
    if (periodic_axis[0] == false && periodic_axis[1] == false) {
        shift_image_domain.push_back(F64vec(0.0));  // NOTE: sign is plus
        return;
    }
#else
    if (periodic_axis[0] == false && periodic_axis[1] == false && periodic_axis[2] == false) {
        shift_image_domain.push_back(F64vec(0.0));  // NOTE: sign is plus
        return;
    }
#endif
    S32ort pos_my_domain_idx, pos_target_domain_idx;
    pos_my_domain_idx.low_ = morton_key.getPMCellID(pos_my_domain.low_);
    pos_my_domain_idx.high_ = morton_key.getPMCellID(pos_my_domain.high_);
    pos_target_domain_idx.low_ = morton_key.getPMCellID(pos_target_domain.low_) - S32vec(icut);
    pos_target_domain_idx.high_ = morton_key.getPMCellID(pos_target_domain.high_) + S32vec(icut);
    if (pos_my_domain_idx.overlapped(pos_target_domain_idx)) {
        shift_image_domain.push_back(F64vec(0.0));  // NOTE: sign is plus
    }
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
    for (S32 lev = 1;; lev++) {
        S32 n_image_per_level = 0;
        S32 lev_x = 0;
        S32 lev_y = 0;
        if (periodic_axis[0] == true) lev_x = lev;
        if (periodic_axis[1] == true) lev_y = lev;
        for (S32 ix = -lev_x; ix <= lev_x; ix++) {
            for (S32 iy = -lev_y; iy <= lev_y; iy++) {
                if ((std::abs(ix) != lev_x || periodic_axis[0] == false) && (std::abs(iy) != lev_y || periodic_axis[1] == false)) continue;
                const F64vec shift_tmp(ix * size_root_domain.x, iy * size_root_domain.y);
                const F64ort pos_image_tmp = pos_target_domain.shift(shift_tmp);
                S32ort pos_image_tmp_idx;
                pos_image_tmp_idx.low_ = morton_key.getPMCellID(pos_image_tmp.low_) - S32vec(icut);
                pos_image_tmp_idx.high_ = morton_key.getPMCellID(pos_image_tmp.high_) + S32vec(icut);
                if (pos_my_domain_idx.overlapped(pos_image_tmp_idx)) {
                    shift_image_domain.push_back(shift_tmp);  // NOTE: sign is plus
                    n_image_per_level++;
                }
            }
        }
        if (n_image_per_level == 0) break;
    }
#else
    for (S32 lev = 1;; lev++) {
        S32 n_image_per_level = 0;
        S32 lev_x = 0;
        S32 lev_y = 0;
        S32 lev_z = 0;
        if (periodic_axis[0] == true) lev_x = lev;
        if (periodic_axis[1] == true) lev_y = lev;
        if (periodic_axis[2] == true) lev_z = lev;
        for (S32 ix = -lev_x; ix <= lev_x; ix++) {
            for (S32 iy = -lev_y; iy <= lev_y; iy++) {
                for (S32 iz = -lev_z; iz <= lev_z; iz++) {
                    if ((std::abs(ix) != lev_x || periodic_axis[0] == false) && (std::abs(iy) != lev_y || periodic_axis[1] == false) &&
                        (std::abs(iz) != lev_z || periodic_axis[2] == false))
                        continue;
                    const F64vec shift_tmp(ix * size_root_domain.x, iy * size_root_domain.y, iz * size_root_domain.z);
                    const F64ort pos_image_tmp = pos_target_domain.shift(shift_tmp);
                    S32ort pos_image_tmp_idx;
                    pos_image_tmp_idx.low_ = morton_key.getPMCellID(pos_image_tmp.low_) - S32vec(icut);
                    pos_image_tmp_idx.high_ = morton_key.getPMCellID(pos_image_tmp.high_) + S32vec(icut);
                    if (pos_my_domain_idx.overlapped(pos_image_tmp_idx)) {
                        shift_image_domain.push_back(shift_tmp);  // NOTE: sign is plus
                        n_image_per_level++;
                    }
                }
            }
        }
        if (n_image_per_level == 0) break;
    }
#endif
}

static inline S32vec CalcIDOfImageDomain(const F64ort &pos_root_domain, const F64vec &pos_target, const bool periodic_axis[]) {
    if (pos_root_domain.contained(pos_target)) {
        return S32vec(0.0);
    }
    const F64vec size_root_domain = pos_root_domain.getFullLength();
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
    for (S32 lev = 1;; lev++) {
        S32 lev_x = 0;
        S32 lev_y = 0;
        if (periodic_axis[0] == true) lev_x = lev;
        if (periodic_axis[1] == true) lev_y = lev;
        for (S32 ix = -lev_x; ix <= lev_x; ix++) {
            for (S32 iy = -lev_y; iy <= lev_y; iy++) {
                if ((std::abs(ix) != lev_x || periodic_axis[0] == false) && (std::abs(iy) != lev_y || periodic_axis[1] == false)) continue;
                const F64vec shift_tmp(ix * size_root_domain.x, iy * size_root_domain.y);
                if (pos_root_domain.contained(pos_target + shift_tmp)) {
                    return S32vec(-ix, -iy);
                }
            }
        }
    }
#else
    for (S32 lev = 1;; lev++) {
        S32 lev_x = 0;
        S32 lev_y = 0;
        S32 lev_z = 0;
        if (periodic_axis[0] == true) lev_x = lev;
        if (periodic_axis[1] == true) lev_y = lev;
        if (periodic_axis[2] == true) lev_z = lev;
        for (S32 ix = -lev_x; ix <= lev_x; ix++) {
            for (S32 iy = -lev_y; iy <= lev_y; iy++) {
                for (S32 iz = -lev_z; iz <= lev_z; iz++) {
                    if ((std::abs(ix) != lev_x || periodic_axis[0] == false) && (std::abs(iy) != lev_y || periodic_axis[1] == false) &&
                        (std::abs(iz) != lev_z || periodic_axis[2] == false))
                        continue;
                    const F64vec shift_tmp(ix * size_root_domain.x, iy * size_root_domain.y, iz * size_root_domain.z);
                    if (pos_root_domain.contained(pos_target + shift_tmp)) {
                        return S32vec(-ix, -iy, -iz);
                    }
                }
            }
        }
    }
#endif
}

// it only works for DIMENSION = 3 or 2
static PS_INLINE S32vec GetRankVec(const S32 n_domain[DIMENSION_LIMIT], const S32 rank_glb) {
    S32vec rank_1d;
    S32 rank_tmp = rank_glb;
    S32 div = 1;
    for (auto i = DIMENSION - 1; i >= 0; i--) {
        rank_1d[i] = rank_tmp % n_domain[i];
        div *= n_domain[i];
        rank_tmp = rank_glb / div;
    }
    return rank_1d;
}

// it also works for arbitrary dimension
static PS_INLINE void GetRankVec(S32 *rank_vec, const S32 *n_domain, const S32 rank_glb, const S32 dim) {
    S32 rank_tmp = rank_glb;
    S32 div = 1;
    for (auto i = dim - 1; i >= 0; i--) {
        rank_vec[i] = rank_tmp % n_domain[i];
        div *= n_domain[i];
        rank_tmp = rank_glb / div;
    }
}

// it works for dim = 3 or 2
static PS_INLINE S32 GetRankGlb(const S32vec rank_vec, const S32 n_domain[DIMENSION_LIMIT]) {
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    return rank_vec[0] * n_domain[1] + rank_vec[1];
#else
    return rank_vec[0] * n_domain[1] * n_domain[2] + rank_vec[1] * n_domain[2] + rank_vec[2];
#endif
}
// it works for arbitrary dimension
static PS_INLINE S32 GetRankGlb(const S32 *rank_vec, const S32 *n_domain, const S32 dim) {
    S32 rank_new = 0;
    S32 radix = 1;
    for (auto i = dim - 1; i >= 0; i--) {
        rank_new += radix * rank_vec[i];
        radix *= n_domain[i];
    }
    return rank_new;
}

// it works for dim = 3 or 2
static PS_INLINE S32 GetRank1d(const S32 rank_glb, const S32 n_domain[DIMENSION_LIMIT], const S32 cid) {
    S32 ret = -1;
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    if (cid == 0)
        ret = rank_glb / n_domain[1];
    else if (cid == 1)
        ret = rank_glb % n_domain[1];
#else
    if (cid == 0)
        ret = (rank_glb / n_domain[2]) / n_domain[1];
    else if (cid == 1)
        ret = (rank_glb / n_domain[2]) % n_domain[1];
    else if (cid == 2)
        ret = rank_glb % n_domain[2];
#endif
    return ret;
}
// it works for arbitrary dimension
static PS_INLINE S32 GetRank1d(const S32 rank_glb, const S32 *n_domain, const S32 dim, const S32 cid) {
    S32 rank_new = rank_glb;
    for (auto i = dim - 1; i > cid; i--) {
        rank_new /= n_domain[i];
    }
    return rank_new % n_domain[cid];
}

static PS_INLINE S32 GetAppropriateRank1D(const F64vec &pos, const S32 n_domain[DIMENSION_LIMIT], const F64ort domain[], const int cid,
                                          const S32 rank_glb_search_start = 0) {
    const auto pos_1d = pos[cid];
    const auto rank_vec_org = GetRankVec(n_domain, rank_glb_search_start);
    const auto rank_1d_org = rank_vec_org[cid];
    const auto rank_1d_mid = (rank_1d_org + n_domain[cid] + n_domain[cid] / 2) % n_domain[cid];
    auto rank_vec_mid = rank_vec_org;
    rank_vec_mid[cid] = rank_1d_mid;
    const auto rank_glb_mid = GetRankGlb(rank_vec_mid, n_domain);
    S32 shift = 1;
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    if (cid == 0) shift = n_domain[1];
#else
    if (cid == 1)
        shift = n_domain[2];
    else if (cid == 0)
        shift = n_domain[1] * n_domain[2];
#endif
    auto rank_vec_start = rank_vec_org;
    if (rank_1d_org > rank_1d_mid && domain[rank_glb_mid].low_[cid] > pos_1d) {
        rank_vec_start[cid] = 0;
    } else if (rank_1d_org < rank_1d_mid && domain[rank_glb_mid].high_[cid] <= pos_1d) {
        rank_vec_start[cid] = n_domain[cid] - 1;
    }
    auto rank_glb_tmp = GetRankGlb(rank_vec_start, n_domain);
    // int n_shift = 0;
    while (pos_1d >= domain[rank_glb_tmp].high_[cid]) {
        rank_glb_tmp += shift;
        // n_shift++;
    }
    while (pos_1d < domain[rank_glb_tmp].low_[cid]) {
        rank_glb_tmp -= shift;
        // n_shift++;
    }
    // std::cerr<<"n_shift= "<<n_shift<<std::endl;
    return GetRank1d(rank_glb_tmp, n_domain, cid);
}

inline void CalcNumShift(const F64ort root_domain, const F64ort my_outer_boundary, const F64vec shift, S32 &n_shift) {
    F64ort root_domain_tmp = root_domain;
    root_domain_tmp.high_ += shift;
    root_domain_tmp.low_ += shift;
    while (my_outer_boundary.overlapped(root_domain_tmp)) {
        root_domain_tmp.high_ += shift;
        root_domain_tmp.low_ += shift;
        n_shift++;
    }
}
}  // namespace ParticleSimulator
#include <utils_comm.hpp>
