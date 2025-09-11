#pragma once

#include<limits>
#include"vector3.hpp"

// T must be S32 or S64.
namespace ParticleSimulator{
    template<class T>
    class Orthotope3i{
    public:
        Vector3<T> low_;
        Vector3<T> high_;

        Orthotope3i()
            : low_(std::numeric_limits<T>::max()), 
              high_(-std::numeric_limits<T>::max()) {}

        Orthotope3i(const Vector3<T> & _low, const Vector3<T> & _high)
            : low_(_low), high_(_high) {}

        Orthotope3i(const Orthotope3i & src) : low_(src.low_), high_(src.high_){}

        const Orthotope3i & operator = (const Orthotope3i & rhs){
            high_ = rhs.high_;
            low_ = rhs.low_;
            return (*this);
        }

        void init(){
            low_  = std::numeric_limits<T>::max();
            high_ = -low_;
        }

        unsigned int isNotValid() const {
            return (this->high_.x < this->low_.x)
                || (this->high_.y < this->low_.y)
                || (this->high_.z < this->low_.z);
        }

        unsigned int isValid() const {
            return isNotValid() ^ 0x1;
        }

        Orthotope3i shift(const Vector3<T> & vec) const {
            return Orthotope3i(low_ + vec, high_ + vec);
        }

        void merge(const Orthotope3i & ort){
            this->high_.x = ( this->high_.x > ort.high_.x ) ? this->high_.x : ort.high_.x;
            this->high_.y = ( this->high_.y > ort.high_.y ) ? this->high_.y : ort.high_.y;
            this->high_.z = ( this->high_.z > ort.high_.z ) ? this->high_.z : ort.high_.z;
            this->low_.x = ( this->low_.x <= ort.low_.x ) ? this->low_.x : ort.low_.x;
            this->low_.y = ( this->low_.y <= ort.low_.y ) ? this->low_.y : ort.low_.y;
            this->low_.z = ( this->low_.z <= ort.low_.z ) ? this->low_.z : ort.low_.z;
        }

        void merge(const Vector3<T> & vec){
            this->high_.x = ( this->high_.x > vec.x ) ? this->high_.x : vec.x;
            this->high_.y = ( this->high_.y > vec.y ) ? this->high_.y : vec.y;
            this->high_.z = ( this->high_.z > vec.z ) ? this->high_.z : vec.z;
            this->low_.x = ( this->low_.x <= vec.x ) ? this->low_.x : vec.x;
            this->low_.y = ( this->low_.y <= vec.y ) ? this->low_.y : vec.y;
            this->low_.z = ( this->low_.z <= vec.z ) ? this->low_.z : vec.z;
        }

        void merge(const Vector3<T> & vec, const T size){
            this->high_.x = ( this->high_.x > vec.x + size ) ? this->high_.x : vec.x + size;
            this->high_.y = ( this->high_.y > vec.y + size ) ? this->high_.y : vec.y + size;
            this->high_.z = ( this->high_.z > vec.z + size) ? this->high_.z : vec.z + size;
            this->low_.x = ( this->low_.x <= vec.x - size) ? this->low_.x : vec.x - size;
            this->low_.y = ( this->low_.y <= vec.y - size) ? this->low_.y : vec.y - size;
            this->low_.z = ( this->low_.z <= vec.z - size) ? this->low_.z : vec.z - size;
        }

        unsigned int isNotContainedBy(const Orthotope3i & a) const {
            const Vector3<T> a_high = a.high_;
            const Vector3<T> a_low = a.low_;
            const Vector3<T> b_high = this->high_;
            const Vector3<T> b_low = this->low_;
            return (b_low.x < a_low.x) || (a_high.x < b_high.x)
                || (b_low.y < a_low.y) || (a_high.y < b_high.y)
                || (b_low.z < a_low.z) || (a_high.z < b_high.z);
        }

        unsigned int isContainedBy(const Orthotope3i & a) const {
            if (a.isValid()) return isNotContainedBy(a) ^ 0x1;
            else return 0;
            // The result is the same as that of the following code:
            //return (a_low.x <= b_low.x) && (b_high.x <= a_high.x)
            //    && (a_low.y <= b_low.y) && (b_high.y <= a_high.y)
            //    && (a_low.z <= b_low.z) && (b_high.z <= a_high.z);
        }

        unsigned int doesNotContain(const Vector3<T> & pos) const {
            return (pos.x < this->low_.x) || (this->high_.x < pos.x)
                || (pos.y < this->low_.y) || (this->high_.y < pos.y)
                || (pos.z < this->low_.z) || (this->high_.z < pos.z);
        }

        unsigned int contains(const Vector3<T> & pos) const {
            return doesNotContain(pos) ^ 0x1;
            // The result is the same as that of the following code:
            //return (this->low_.x <= pos.x) && (pos.x <= this->high_.x)
            //    && (this->low_.y <= pos.y) && (pos.y <= this->high_.y)
            //    && (this->low_.z <= pos.z) && (pos.z <= this->high_.z);
            // Note that it differs from orthotope3.hpp.
        }

        unsigned int doesNotContain(const Orthotope3i & a) const {
            return (a.low_.x < this->low_.x) || (this->high_.x < a.high_.x)
                || (a.low_.y < this->low_.y) || (this->high_.y < a.high_.y)
                || (a.low_.z < this->low_.z) || (this->high_.z < a.high_.z);
        }

        unsigned int contains(const Orthotope3i & a) const {
            if (a.isValid()) return doesNotContain(a) ^ 0x1;
            else return 0;
            // The result is the same as that of the following code:
            //return (this->low_.x <= a.low_.x) && (a.high_.x <= this->high_.x)
            //    && (this->low_.y <= a.low_.y) && (a.high_.y <= this->high_.y)
            //    && (this->low_.z <= a.low_.z) && (a.high_.z <= this->high_.z);
        }

        unsigned int notOverlapped(const Vector3<T> & pos) const {
            return (pos.x < low_.x) || (high_.x < pos.x)
                || (pos.y < low_.y) || (high_.y < pos.y)
                || (pos.z < low_.z) || (high_.z < pos.z);
        }
        unsigned int overlapped(const Vector3<T> & pos) const {
            return notOverlapped(pos) ^ 0x1;
        }

        unsigned int notOverlapped(const Orthotope3i & a) const {
            return (a.high_.x < low_.x) || (high_.x < a.low_.x)
                || (a.high_.y < low_.y) || (high_.y < a.low_.y)
                || (a.high_.z < low_.z) || (high_.z < a.low_.z);
        }

        unsigned int overlapped(const Orthotope3i & a) const {
            return notOverlapped(a) ^ 0x1;
        }

        const Vector3<T> getLength() const {
            return (high_ - low_ + Vector3<T>(1));
        }

        const Vector3<T> getDistanceEach(const Vector3<T> & vec) const {
            const Vector3<T> b_high = this->high_;
            const Vector3<T> b_low = this->low_;
            T dx = (vec.x > b_high.x) ? (vec.x - b_high.x) : ( (vec.x < b_low.x) ? (b_low.x - vec.x) : T(0) );
            T dy = (vec.y > b_high.y) ? (vec.y - b_high.y) : ( (vec.y < b_low.y) ? (b_low.y - vec.y) : T(0) );
            T dz = (vec.z > b_high.z) ? (vec.z - b_high.z) : ( (vec.z < b_low.z) ? (b_low.z - vec.z) : T(0) );
            return Vector3<T>(dx, dy, dz);
        }

        const T getVolume() const {
            const Vector3<T> len = getLength();
            return len.x * len.y * len.z;
        }
        
        bool calcIntersection(const Orthotope3i & tgt,
                              Orthotope3i & ret) const {
            if (overlapped(tgt)) {
                ret.low_.x  = std::max(this->low_.x, tgt.low_.x);
                ret.high_.x = std::min(this->high_.x, tgt.high_.x);
                ret.low_.y  = std::max(this->low_.y, tgt.low_.y);
                ret.high_.y = std::min(this->high_.y, tgt.high_.y);
                ret.low_.z  = std::max(this->low_.z, tgt.low_.z);
                ret.high_.z = std::min(this->high_.z, tgt.high_.z);
                return true;
            } else {
                return false;
            }
        }

        Orthotope3i getIntersectionNoCheck(const Orthotope3i & tgt) const {
            Orthotope3i ret;
            ret.low_.x  = std::max(this->low_.x, tgt.low_.x);
            ret.high_.x = std::min(this->high_.x, tgt.high_.x);
            ret.low_.y  = std::max(this->low_.y, tgt.low_.y);
            ret.high_.y = std::min(this->high_.y, tgt.high_.y);
            ret.low_.z  = std::max(this->low_.z, tgt.low_.z);
            ret.high_.z = std::min(this->high_.z, tgt.high_.z);
            return ret;
        }
        
        //Cast to Orthotope3<U>
        template <typename U>
        operator Orthotope3i<U> () const {
            return Orthotope3i<U> (static_cast < Vector3<U> > (low_),
                                   static_cast < Vector3<U> > (high_));
        }

        friend std::ostream & operator <<(std::ostream & c, const Orthotope3i & u){
            c<<std::setprecision(15)<<u.low_<<"   "<<u.high_;
            return c;
        }
    };
}

