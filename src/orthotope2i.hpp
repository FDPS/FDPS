#pragma once

#include<limits>
#include"vector3.hpp"

// T must be S32 or S64.
namespace ParticleSimulator{
    template<class T>
    class Orthotope2i{
    public:
        Vector2<T> low_;
        Vector2<T> high_;

        Orthotope2i()
            : low_(std::numeric_limits<T>::max()), 
              high_(-std::numeric_limits<T>::max()) {}

        Orthotope2i(const Vector2<T> & _low, const Vector2<T> & _high)
            : low_(_low), high_(_high) {}

        Orthotope2i(const Orthotope2i & src) : low_(src.low_), high_(src.high_){}

        const Orthotope2i & operator = (const Orthotope2i & rhs){
            high_ = rhs.high_;
            low_ = rhs.low_;
            return (*this);
        }

        void init(){
            low_  = std::numeric_limits<T>::max();
            high_ = -low_;
        }

        Orthotope2i shift(const Vector2<T> & vec) const {
            return Orthotope2i(low_ + vec, high_ + vec);
        }

        void merge(const Orthotope2i & ort){
            this->high_.x = ( this->high_.x > ort.high_.x ) ? this->high_.x : ort.high_.x;
            this->high_.y = ( this->high_.y > ort.high_.y ) ? this->high_.y : ort.high_.y;
            this->low_.x = ( this->low_.x <= ort.low_.x ) ? this->low_.x : ort.low_.x;
            this->low_.y = ( this->low_.y <= ort.low_.y ) ? this->low_.y : ort.low_.y;
        }

        void merge(const Vector2<T> & vec){
            this->high_.x = ( this->high_.x > vec.x ) ? this->high_.x : vec.x;
            this->high_.y = ( this->high_.y > vec.y ) ? this->high_.y : vec.y;
            this->low_.x = ( this->low_.x <= vec.x ) ? this->low_.x : vec.x;
            this->low_.y = ( this->low_.y <= vec.y ) ? this->low_.y : vec.y;
        }

        void merge(const Vector2<T> & vec, const T size){
            this->high_.x = ( this->high_.x > vec.x + size ) ? this->high_.x : vec.x + size;
            this->high_.y = ( this->high_.y > vec.y + size ) ? this->high_.y : vec.y + size;
            this->low_.x = ( this->low_.x <= vec.x - size) ? this->low_.x : vec.x - size;
            this->low_.y = ( this->low_.y <= vec.y - size) ? this->low_.y : vec.y - size;
        }

        unsigned int notOverlapped(const Vector2<T> & pos) const {
            return (pos.x < low_.x) || (high_.x < pos.x)
                || (pos.y < low_.y) || (high_.y < pos.y);
        }
        unsigned int overlapped(const Vector2<T> & pos) const {
            return notOverlapped(pos) ^ 0x1;
        }

        unsigned int notOverlapped(const Orthotope2i & a) const {
            return (a.high_.x < low_.x) || (high_.x < a.low_.x)
                || (a.high_.y < low_.y) || (high_.y < a.low_.y);
        }

        unsigned int overlapped(const Orthotope2i & a) const {
            return notOverlapped(a) ^ 0x1;
        }

        const Vector2<T> getLength() const {
            return (high_ - low_);
        }

        const Vector2<T> getDistanceEach(const Vector2<T> & vec) const {
            const Vector2<T> b_high = this->high_;
            const Vector2<T> b_low = this->low_;
            T dx = (vec.x > b_high.x) ? (vec.x - b_high.x) : ( (vec.x < b_low.x) ? (b_low.x - vec.x) : T(0) );
            T dy = (vec.y > b_high.y) ? (vec.y - b_high.y) : ( (vec.y < b_low.y) ? (b_low.y - vec.y) : T(0) );
            return Vector2<T>(dx, dy);
        }
        
        //Cast to Orthotope2i<U>
        template <typename U>
        operator Orthotope2i<U> () const {
            return Orthotope2i<U> (static_cast < Vector2<U> > (low_),
                                   static_cast < Vector2<U> > (high_));
        }

        friend std::ostream & operator <<(std::ostream & c, const Orthotope2i & u){
            c<<std::setprecision(15)<<u.low_<<"   "<<u.high_;
            return c;
        }
    };
}
