#pragma once

#include<iostream>
#include<iomanip>
#include<cmath>
#include"vector2.hpp"

namespace ParticleSimulator{
    template<class T>
    class Matrix2{
    public:
        //constructor
        T xx, yy, xy, yx;
        Matrix2() : xx(T(0)), yy(T(0)), xy(T(0)), yx(T(0)) {} 
        Matrix2(const T _xx, const T _yy, const T _xy, const T _yx) 
            : xx(_xx), yy(_yy), xy(_xy), yx(_yx) {} 
        Matrix2(const T s) : xx(s), yy(s), xy(s), yx(s){}
        Matrix2(const Matrix2 & src) : xx(src.xx), yy(src.yy), xy(src.xy), yx(src.yx) {}
        Matrix2(const Vector2<T> & a, const Vector2<T> & b) : xx(a.x * b.x), yy(a.y * b.y), xy(a.x * b.y), yx(a.y * b.x) {}

        const Matrix2 & operator = (const Matrix2 & rhs) {
            xx = rhs.xx;
            yy = rhs.yy;
            xy = rhs.xy;
            yx = rhs.yx;
            return (*this);
        }
        const Matrix2 & operator = (const T s) {
            xx = yy = xy = yx = s;
            return (*this);
        }

        Matrix2 operator + (const Matrix2 & rhs) const {
            return Matrix2(xx + rhs.xx, yy + rhs.yy, xy + rhs.xy, yx + rhs.yx);
        }
        Matrix2 operator + (const T & rhs) const {
            return Matrix2(xx + rhs, yy + rhs, xy, yx);
        }
        Matrix2 operator - (const T & rhs) const {
            return Matrix2(xx - rhs, yy - rhs, xy, yx);
        }
        const Matrix2 & operator += (const Matrix2 & rhs) {
            (*this) = (*this) + rhs;
            return (*this);
        }
        Matrix2 operator - (const Matrix2 & rhs) const {
            return Matrix2(xx - rhs.xx, yy - rhs.yy, xy - rhs.xy, yx - rhs.yx);
        }
        const Matrix2 & operator -= (const Matrix2 & rhs) {
            (*this) = (*this) - rhs;
            return (*this);
        }
        Matrix2 operator * (const T & rhs) const {
            return Matrix2(xx * rhs, yy * rhs, xy * rhs, yx * rhs);
        }
        Vector2<T> operator * (const Vector2<T> & rhs) const {
            return Vector2<T>(xx * rhs.x + xy * rhs.y, yx * rhs.x + yy * rhs.y);
        }
        Matrix2 operator * (const Matrix2 & rhs) const {
            return Matrix2(xx * rhs.xx + xy * rhs.yx, yx * rhs.xy + yy * rhs.yy, xx + rhs.xy + xy * rhs.yy, yx * rhs.xx + yy * rhs.yx);
        }
        Matrix2 operator / (const T & rhs) const {
            return Matrix2(xx / rhs, yy / rhs, xy / rhs, yx / rhs);
        }
        const Matrix2 & operator *= (const T & rhs) {
            (*this) = (*this) * rhs;
            return (*this);
        }
        friend Matrix2 operator * (const T s, const Matrix2 & m) {
            return (m * s);
        }
        const Matrix2 & operator /= (const T & rhs) {
            (*this).xx /= rhs;
            (*this).yy /= rhs;
            (*this).xy /= rhs;
            (*this).yx /= rhs;
            return (*this);
        }

        const Matrix2 & operator + () const {
            return (* this);
        }

        const Matrix2 operator - () const {
            return Matrix2(-xx, -yy, -xy, -yx);
        }


        T getTrace() const {
            return (xx + yy);
        }
        T getDeterminant() const {
            return (xx * yy - xy * yx);
        }
        T getSecondInvariant() const {
            return 0.5 * (- this->getTrace() * this->getTrace() + (*this * *this).getTrace());
        }


        template <typename U>
        operator Matrix2<U> () const {
            return Matrix2<U>( static_cast<U>(xx), static_cast<U>(yy), static_cast<U>(xy), static_cast<U>(yx) );
        }

        friend std::ostream& operator << (std::ostream& c, const Matrix2<T> & mat){
            c<<mat.xx<<"   "<<mat.xy<<std::endl;
            c<<mat.yx<<"   "<<mat.yy<<std::endl;
            return c;
        }
    };
}
