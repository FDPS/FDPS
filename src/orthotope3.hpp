#pragma once

#include<limits>
#include"vector3.hpp"

// T must be F32 or F64.
namespace ParticleSimulator{
    template<class T>
    class Orthotope3{
    public:
        Vector3<T> low_;
        Vector3<T> high_;

        Orthotope3(): low_(9999.9), high_(-9999.9){}

        Orthotope3(const Vector3<T> & _low, const Vector3<T> & _high)
            : low_(_low), high_(_high) {}

        Orthotope3(const Orthotope3 & src) : low_(src.low_), high_(src.high_){}

        Orthotope3(const Vector3<T> & center, const T length) :
            low_(center-(Vector3<T>)(length)), high_(center+(Vector3<T>)(length)) {
        }

        const Orthotope3 & operator = (const Orthotope3 & rhs){
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

        Orthotope3 shift( const Vector3<T> & vec) const {
            return Orthotope3(low_ + vec, high_ + vec);
        }

        void merge( const Orthotope3 & ort ){
            this->high_.x = ( this->high_.x > ort.high_.x ) ? this->high_.x : ort.high_.x;
            this->high_.y = ( this->high_.y > ort.high_.y ) ? this->high_.y : ort.high_.y;
            this->high_.z = ( this->high_.z > ort.high_.z ) ? this->high_.z : ort.high_.z;
            this->low_.x = ( this->low_.x <= ort.low_.x ) ? this->low_.x : ort.low_.x;
            this->low_.y = ( this->low_.y <= ort.low_.y ) ? this->low_.y : ort.low_.y;
            this->low_.z = ( this->low_.z <= ort.low_.z ) ? this->low_.z : ort.low_.z;
        }

        void merge( const Vector3<T> & vec ){
            this->high_.x = ( this->high_.x > vec.x ) ? this->high_.x : vec.x;
            this->high_.y = ( this->high_.y > vec.y ) ? this->high_.y : vec.y;
            this->high_.z = ( this->high_.z > vec.z ) ? this->high_.z : vec.z;
            this->low_.x = ( this->low_.x <= vec.x ) ? this->low_.x : vec.x;
            this->low_.y = ( this->low_.y <= vec.y ) ? this->low_.y : vec.y;
            this->low_.z = ( this->low_.z <= vec.z ) ? this->low_.z : vec.z;
        }

        void merge( const Vector3<T> & vec, const T size){
            this->high_.x = ( this->high_.x > vec.x + size ) ? this->high_.x : vec.x + size;
            this->high_.y = ( this->high_.y > vec.y + size ) ? this->high_.y : vec.y + size;
            this->high_.z = ( this->high_.z > vec.z + size) ? this->high_.z : vec.z + size;
            this->low_.x = ( this->low_.x <= vec.x - size) ? this->low_.x : vec.x - size;
            this->low_.y = ( this->low_.y <= vec.y - size) ? this->low_.y : vec.y - size;
            this->low_.z = ( this->low_.z <= vec.z - size) ? this->low_.z : vec.z - size;
        }

        unsigned int notContained(const Vector3<T> & pos) const {
            return (pos.x < low_.x) || (high_.x <= pos.x)
                || (pos.y < low_.y) || (high_.y <= pos.y)
                || (pos.z < low_.z) || (high_.z <= pos.z);
        }

        unsigned int contained(const Vector3<T> & pos) const {
            return notContained(pos) ^ 0x1;
        }

        unsigned int notContained(const Orthotope3 & a) const {
            const Vector3<T> a_high = a.high_;
            const Vector3<T> a_low = a.low_;
            const Vector3<T> b_high = this->high_;
            const Vector3<T> b_low = this->low_;

            return (a_high.x < b_low.x) || (b_high.x <= a_low.x)
                || (a_high.y < b_low.y) || (b_high.y <= a_low.y)
                || (a_high.z < b_low.z) || (b_high.z <= a_low.z);
        }

        unsigned int contained(const Orthotope3 & a) const {
            return notContained(a) ^ 0x1;
        }

        unsigned int notOverlapped(const Vector3<T> & pos) const {
            return (pos.x < low_.x) || (high_.x < pos.x)
                || (pos.y < low_.y) || (high_.y < pos.y)
                || (pos.z < low_.z) || (high_.z < pos.z);
        }
        unsigned int overlapped(const Vector3<T> & pos) const {
            return notOverlapped(pos) ^ 0x1;
        }

        unsigned int notOverlapped(const Orthotope3 & a) const {
            return (a.high_.x < low_.x) || (high_.x < a.low_.x)
                || (a.high_.y < low_.y) || (high_.y < a.low_.y)
                || (a.high_.z < low_.z) || (high_.z < a.low_.z);
        }

        unsigned int overlapped(const Orthotope3 & a) const {
            return notOverlapped(a) ^ 0x1;
        }

        const Vector3<T> getCenter() const {
            return (high_ + low_) * T(0.5);
        }
        const Vector3<T> getHalfLength() const {
            return (high_ - low_) * T(0.5);
        }
        const Vector3<T> getFullLength() const {
            return (high_ - low_);
        }

        const T getDistanceMinSQ(const Orthotope3 & ort) const {
            // it dose not work in the case of an orthotope with negative volume
            const Vector3<T> a_high = ort.high_;
            const Vector3<T> a_low = ort.low_;
            const Vector3<T> b_high = this->high_;
            const Vector3<T> b_low = this->low_;
            const T x0 = b_low.x - a_high.x;
            const T x1 = a_low.x - b_high.x;
            const T y0 = b_low.y - a_high.y;
            const T y1 = a_low.y - b_high.y;
            const T z0 = b_low.z - a_high.z;
            const T z1 = a_low.z - b_high.z;
            T dx = (x0 > x1) ? x0 : x1;
            T dy = (y0 > y1) ? y0 : y1;
            T dz = (z0 > z1) ? z0 : z1;
            dx = (x0*x1 < 0) ? dx : T(0);
            dy = (y0*y1 < 0) ? dy : T(0);
            dz = (z0*z1 < 0) ? dz : T(0);
            return dx*dx + dy*dy + dz*dz;
        }
        const T getDistanceMinSq(const Orthotope3 & ort) const {
            return getDistanceMinSQ(ort);
        }
        
        const T getDistanceMinSQ(const Vector3<T> & vec) const {
            const Vector3<T> b_high = this->high_;
            const Vector3<T> b_low = this->low_;
            T dx = (vec.x > b_high.x) ? (vec.x - b_high.x) : ( (vec.x < b_low.x) ? (b_low.x - vec.x) : T(0) );
            T dy = (vec.y > b_high.y) ? (vec.y - b_high.y) : ( (vec.y < b_low.y) ? (b_low.y - vec.y) : T(0) );
            T dz = (vec.z > b_high.z) ? (vec.z - b_high.z) : ( (vec.z < b_low.z) ? (b_low.z - vec.z) : T(0) );
            return dx*dx + dy*dy + dz*dz;
        }

        const T getDistanceMinSq(const Vector3<T> & vec) const {
            return getDistanceMinSQ(vec);
        }
        
        //Cast to Orthotope3<U>
        template <typename U>
        operator Orthotope3<U> () const {
            return Orthotope3<U> (static_cast < Vector3<U> > (low_),
                                  static_cast < Vector3<U> > (high_));
        }

        friend std::ostream & operator <<(std::ostream & c, const Orthotope3 & u){
            c<<std::setprecision(15)<<u.low_<<"   "<<u.high_;
            return c;
        }
    };
}
