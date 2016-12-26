#pragma once

#include<iostream>
#include<iomanip>
#include<cmath>
#include"vector3.hpp"

namespace ParticleSimulator{
    template<class T>
    class MatrixSym3{
    public:
        //constructor
        T xx, yy, zz, xy, xz, yz;
        MatrixSym3() : xx(T(0)), yy(T(0)), zz(T(0)), xy(T(0)), xz(T(0)), yz(T(0)) {} 
        MatrixSym3(const T _xx, const T _yy, const T _zz, 
                   const T _xy, const T _xz, const T _yz ) 
            : xx(_xx), yy(_yy), zz(_zz), xy(_xy), xz(_xz), yz(_yz) {} 
        MatrixSym3(const T s) : xx(s), yy(s), zz(s), xy(s), xz(s), yz(s) {}
        MatrixSym3(const MatrixSym3 & src) : xx(src.xx), yy(src.yy), zz(src.zz), 
                                             xy(src.xy), xz(src.xz), yz(src.yz) {}

        const MatrixSym3 & operator = (const MatrixSym3 & rhs) {
            xx = rhs.xx;
            yy = rhs.yy;
            zz = rhs.zz;
            xy = rhs.xy;
            xz = rhs.xz;
            yz = rhs.yz;
            return (*this);
        }
        const MatrixSym3 & operator = (const T s) {
            xx = yy = zz = xy = xz = yz = s;
            return (*this);
        }

        MatrixSym3 operator + (const MatrixSym3 & rhs) const {
            return MatrixSym3( xx + rhs.xx, yy + rhs.yy, zz + rhs.zz,
                               xy + rhs.xy, xz + rhs.xz, yz + rhs.yz);
        }
        const MatrixSym3 & operator += (const MatrixSym3 & rhs) {
            (*this) = (*this) + rhs;
            return (*this);
        }
        MatrixSym3 operator - (const MatrixSym3 & rhs) const {
            return MatrixSym3( xx - rhs.xx, yy - rhs.yy, zz - rhs.zz,
                               xy - rhs.xy, xz - rhs.xz, yz - rhs.yz);
        }
        const MatrixSym3 & operator -= (const MatrixSym3 & rhs) {
            (*this) = (*this) - rhs;
            return (*this);
        }

        T getTrace() const {
            return (xx + yy + zz);
        }

        template <typename U>
        operator MatrixSym3<U> () const {
            return MatrixSym3<U>( static_cast<U>(xx), static_cast<U>(yy), static_cast<U>(zz),
                                  static_cast<U>(xy), static_cast<U>(xz), static_cast<U>(yz) );
        }

        friend std::ostream & operator << (std::ostream & c, const MatrixSym3 & mat){
            c<<mat.xx<<"   "<<mat.xy<<"    "<<mat.xz<<std::endl;
            c<<mat.xy<<"   "<<mat.yy<<"    "<<mat.yz<<std::endl;
            c<<mat.xz<<"   "<<mat.yz<<"    "<<mat.zz<<std::endl;
            return c;
        }
    };
}
