#pragma once

#include<limits>
#include"vector2.hpp"

// T must be F32 or F64.
namespace ParticleSimulator{
    template<class T>
    class Orthotope2{
    public:
        Vector2<T> low_;
        Vector2<T> high_;

        Orthotope2(): low_(9999.9), high_(-9999.9){}
        
        Orthotope2(const Vector2<T> & _low, const Vector2<T> & _high)
            : low_(_low), high_(_high){}
        
        Orthotope2(const Orthotope2 & src) : low_(src.low_), high_(src.high_){}

        Orthotope2(const Vector2<T> & center, const T length) :
            low_(center-(Vector2<T>)(length)), high_(center+(Vector2<T>)(length)) {
        }

        const Orthotope2 & operator = (const Orthotope2 & rhs){
            high_ = rhs.high_;
            low_ = rhs.low_;
            return (*this);
        }

        void initNegativeVolume(){
            low_ = std::numeric_limits<float>::max() / 128;
            high_ = -low_;
        }

        void init(){
            initNegativeVolume();
        }

        unsigned int isNotValid() const { // for PMMM
            return (this->high_.x < this->low_.x)
                || (this->high_.y < this->low_.y);
        }

        unsigned int isValid() const { // for PMMM
            return isNotValid() ^ 0x1;
        }


        Orthotope2 shift( const Vector2<T> & vec) const {
            return Orthotope2(low_ + vec, high_ + vec);
        }

        void merge( const Orthotope2 & ort ){
            this->high_.x = ( this->high_.x > ort.high_.x ) ? this->high_.x : ort.high_.x;
            this->high_.y = ( this->high_.y > ort.high_.y ) ? this->high_.y : ort.high_.y;
            this->low_.x = ( this->low_.x <= ort.low_.x ) ? this->low_.x : ort.low_.x;
            this->low_.y = ( this->low_.y <= ort.low_.y ) ? this->low_.y : ort.low_.y;
        }

        void merge( const Vector2<T> & vec ){
            this->high_.x = ( this->high_.x > vec.x ) ? this->high_.x : vec.x;
            this->high_.y = ( this->high_.y > vec.y ) ? this->high_.y : vec.y;
            this->low_.x = ( this->low_.x <= vec.x ) ? this->low_.x : vec.x;
            this->low_.y = ( this->low_.y <= vec.y ) ? this->low_.y : vec.y;
        }

        void merge( const Vector2<T> & vec, const T size){
            this->high_.x = ( this->high_.x > vec.x + size ) ? this->high_.x : vec.x + size;
            this->high_.y = ( this->high_.y > vec.y + size ) ? this->high_.y : vec.y + size;
            this->low_.x = ( this->low_.x <= vec.x - size) ? this->low_.x : vec.x - size;
            this->low_.y = ( this->low_.y <= vec.y - size) ? this->low_.y : vec.y - size;
        }

#if 0
/*
        unsigned int notOverlapped(const Orthotope2 & a) const {
            return (a.high_.x < low_.x) || (high_.x <= a.low_.x)
                || (a.high_.y < low_.y) || (high_.y <= a.low_.y);
        }
*/

        unsigned int notOverlapped(const Orthotope2 & a) const {
            return (a.high_.x < low_.x) || (high_.x < a.low_.x)
                || (a.high_.y < low_.y) || (high_.y < a.low_.y);
        }

        unsigned int overlapped(const Orthotope2 & a) const {
            return notOverlapped(a) ^ 0x1;
        }

        unsigned int notOverlapped(const Vector2<T> & pos) const {
            return (pos.x < low_.x) || (high_.x <= pos.x)
                || (pos.y < low_.y) || (high_.y <= pos.y);
        }

        unsigned int overlapped(const Vector2<T> & pos) const {
            return notOverlapped(pos) ^ 0x1;
        }

        unsigned int notContains(const Orthotope2 & a) const {
            return (a.low_.x < low_.x) || (high_.x < a.high_.x)
                || (a.low_.y < low_.y) || (high_.y < a.high_.y);
        }

        unsigned int contains(const Orthotope2 & a){
            return notContains(a) ^ 0x1;
        }
#endif

        unsigned int notContained(const Vector2<T> & pos) const {
            return (pos.x < low_.x) || (high_.x <= pos.x)
                || (pos.y < low_.y) || (high_.y <= pos.y);
        }
        unsigned int contained(const Vector2<T> & pos) const {
            return notContained(pos) ^ 0x1;
        }

        unsigned int notContained(const Orthotope2 & a) const {
            const Vector2<T> a_high = a.high_;
            const Vector2<T> a_low = a.low_;
            const Vector2<T> b_high = this->high_;
            const Vector2<T> b_low = this->low_;

            return (a_high.x < b_low.x) || (b_high.x < a_low.x)
                || (a_high.y < b_low.y) || (b_high.y < a_low.y);
        }

        unsigned int contained(const Orthotope2 & a) const {
            return notContained(a) ^ 0x1;
        }

        unsigned int isNotContainedBy(const Orthotope2 & a) const { // for PMMM
            const Vector2<T> a_high = a.high_;
            const Vector2<T> a_low = a.low_;
            const Vector2<T> b_high = this->high_;
            const Vector2<T> b_low = this->low_;
            return (b_low.x < a_low.x) || (a_high.x < b_high.x)
                || (b_low.y < a_low.y) || (a_high.y < b_high.y);
        }

        unsigned int isContainedBy(const Orthotope2 & a) const { // for PMMM
            if (a.isValid()) return isNotContainedBy(a) ^ 0x1;
            return 0;
            // The result is the same as that of the following code:
            //return (a_low.x <= b_low.x) && (b_high.x <= a_high.x)
            //    && (a_low.y <= b_low.y) && (b_high.y <= a_high.y);
        }

        unsigned int doesNotContain(const Vector2<T> & pos) const { // for PMMM
            return (pos.x < this->low_.x) || (this->high_.x <= pos.x)
                || (pos.y < this->low_.y) || (this->high_.y <= pos.y);
        }

        unsigned int contains(const Vector2<T> & pos) const { // for PMMM
            return doesNotContain(pos) ^ 0x1;
            // The result is the same as that of the following code:
            //return (this->low_.x <= pos.x) && (pos.x < this->high_.x)
            //    && (this->low_.y <= pos.y) && (pos.y < this->high_.y);
        }

        unsigned int doesNotContain(const Orthotope2 & a) const { // for PMMM
            return (a.low_.x < this->low_.x) || (this->high_.x < a.high_.x)
                || (a.low_.y < this->low_.y) || (this->high_.y < a.high_.y);
        }

        unsigned int contains(const Orthotope2 & a) const { // for PMMM
            if (a.isValid()) return doesNotContain(a) ^ 0x1;
            else return 0;
            // The result is the same as that of the following code:
            //return (this->low_.x <= a.low_.x) && (a.high_.x <= this->high_.x)
            //    && (this->low_.y <= a.low_.y) && (a.high_.y <= this->high_.y);
        }


        unsigned int notOverlapped(const Vector2<T> & pos) const {
            return (pos.x < low_.x) || (high_.x < pos.x)
                || (pos.y < low_.y) || (high_.y < pos.y);
        }
        unsigned int overlapped(const Vector2<T> & pos) const {
            return notOverlapped(pos) ^ 0x1;
        }

        unsigned int notOverlapped(const Orthotope2 & a) const {
            return (a.high_.x < low_.x) || (high_.x < a.low_.x)
                || (a.high_.y < low_.y) || (high_.y < a.low_.y);
        }

        unsigned int overlapped(const Orthotope2 & a) const {
            return notOverlapped(a) ^ 0x1;
        }

        const Vector2<T> getCenter() const {
            return (high_ + low_) * T(0.5);
        }

        const Vector2<T> getHalfLength() const {
            return (high_ - low_) * T(0.5);
        }

        const Vector2<T> getFullLength() const {
            return (high_ - low_);
        }

        const T getDistanceMinSQ(const Orthotope2 & ort) const {
            T x0 = low_.x - ort.high_.x;
            T x1 = ort.low_.x - high_.x;
            T y0 = low_.y - ort.high_.y;
            T y1 = ort.low_.y - high_.y;
            T dx = (x0 > x1) ? x0 : x1;
            T dy = (y0 > y1) ? y0 : y1;
            dx = (x0*x1 < 0) ? dx : T(0);
            dy = (y0*y1 < 0) ? dy : T(0);
            return dx*dx + dy*dy;
        }

        const T getDistanceMinSq(const Orthotope2 & ort) const {
            return getDistanceMinSQ(ort);
        }

        const T getDistanceMinSQ(const Vector2<T> & vec) const {
            T dx = (vec.x > high_.x) ? (vec.x - high_.x) : ( (vec.x < low_.x) ? (low_.x - vec.x) : T(0) );
            T dy = (vec.y > high_.y) ? (vec.y - high_.y) : ( (vec.y < low_.y) ? (low_.y - vec.y) : T(0) );
            return dx*dx + dy*dy;
        }

        const T getDistanceMinSq(const Vector2<T> & vec) const {
            return getDistanceMinSQ(vec);
        }

        const T getVolume() const {
            const Vector2<T> len = getFullLength();
            return len.x * len.y;
        }

        bool calcVolume(T & ret) const {
            if (isValid()) {
                ret = getVolume();
                return true;
            } else {
                ret = 0.0;
                return false;
            }
        }

        bool calcIntersection(const Orthotope2 & tgt,
                              Orthotope2 & ret) const {
            if (overlapped(tgt)) {
                ret.low_.x  = std::max(this->low_.x, tgt.low_.x);
                ret.high_.x = std::min(this->high_.x, tgt.high_.x);
                ret.low_.y  = std::max(this->low_.y, tgt.low_.y);
                ret.high_.y = std::min(this->high_.y, tgt.high_.y);
                return true;
            } else {
                return false;
            }
        }

        template <typename U>
        operator Orthotope2<U> () const {
            return Orthotope2<U> (static_cast < Vector2<U> > (low_),
                                  static_cast < Vector2<U> > (high_));
        }

        friend std::ostream & operator <<(std::ostream & c, const Orthotope2 & u){
            c<<std::setprecision(15)<<u.low_<<"   "<<u.high_;
            return c;
        }
    };
}
