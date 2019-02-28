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
            return (high_ - low_);
        }

        const Vector3<T> getDistanceEach(const Vector3<T> & vec) const {
            const Vector3<T> b_high = this->high_;
            const Vector3<T> b_low = this->low_;
            T dx = (vec.x > b_high.x) ? (vec.x - b_high.x) : ( (vec.x < b_low.x) ? (b_low.x - vec.x) : T(0) );
            T dy = (vec.y > b_high.y) ? (vec.y - b_high.y) : ( (vec.y < b_low.y) ? (b_low.y - vec.y) : T(0) );
            T dz = (vec.z > b_high.z) ? (vec.z - b_high.z) : ( (vec.z < b_low.z) ? (b_low.z - vec.z) : T(0) );
            return Vector3<T>(dx, dy, dz);
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

