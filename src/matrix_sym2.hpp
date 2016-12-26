#pragma once

#include<iostream>
#include<iomanip>
#include<cmath>
#include"vector2.hpp"

namespace ParticleSimulator{
    template<class T>
    class MatrixSym2{
    public:
        //constructor
        T xx, yy, xy;
        MatrixSym2() : xx(T(0)), yy(T(0)), xy(T(0)) {} 
        MatrixSym2(const T _xx, const T _yy, const T _xy) 
            : xx(_xx), yy(_yy), xy(_xy) {} 
        MatrixSym2(const T s) : xx(s), yy(s), xy(s){}
        MatrixSym2(const MatrixSym2 & src) : xx(src.xx), yy(src.yy), xy(src.xy) {}
        MatrixSym2(const Vector2<T> & a, const Vector2<T> & b) : xx(a.x * b.x), yy(a.y * b.y), xy(a.x * b.y) {}

        const MatrixSym2 & operator = (const MatrixSym2 & rhs) {
            xx = rhs.xx;
            yy = rhs.yy;
            xy = rhs.xy;
            return (*this);
        }
        const MatrixSym2 & operator = (const T s) {
            xx = yy = xy = s;
            return (*this);
        }

        MatrixSym2 operator + (const MatrixSym2 & rhs) const {
            return MatrixSym2(xx + rhs.xx, yy + rhs.yy, xy + rhs.xy);
        }
        const MatrixSym2 & operator += (const MatrixSym2 & rhs) {
            (*this) = (*this) + rhs;
            return (*this);
        }
        MatrixSym2 operator - (const MatrixSym2 & rhs) const {
            return MatrixSym2(xx - rhs.xx, yy - rhs.yy, xy - rhs.xy);
        }
        const MatrixSym2 & operator -= (const MatrixSym2 & rhs) {
            (*this) = (*this) - rhs;
            return (*this);
        }
        MatrixSym2 operator * (const T & rhs) const {
            return MatrixSym2(xx * rhs, yy * rhs, xy * rhs);
        }
        friend MatrixSym2 operator * (const T s, const MatrixSym2 & m) {
            return (m * s);
        }
        const MatrixSym2 & operator /= (const T & rhs) {
            (*this).xx /= rhs;
            (*this).yy /= rhs;
            (*this).xy /= rhs;
            return (*this);
        }


        T getTrace() const {
            return (xx + yy);
        }

        template <typename U>
        operator MatrixSym2<U> () const {
            return MatrixSym2<U>( static_cast<U>(xx), static_cast<U>(yy), static_cast<U>(xy) );
        }

        friend std::ostream& operator << (std::ostream& c, const MatrixSym2<T> & mat){
            c<<mat.xx<<"   "<<mat.xy<<std::endl;
            c<<mat.xy<<"   "<<mat.yy<<std::endl;
            return c;
        }
    };
}
