#pragma once

#include<iostream>
#include<iomanip>
#include<cmath>

namespace ParticleSimulator{
    template<class T>
    class Matrix3{
    public:
      //constructor
      T xx, xy, xz;
      T yx, yy, yz;
      T zx, zy, zz;
      Matrix3()
	: xx(T(0)), xy(T(0)), xz(T(0)),
	  yx(T(0)), yy(T(0)), yz(T(0)),
	  zx(T(0)), zy(T(0)), zz(T(0)) {}

      Matrix3(const T _xx, const T _xy, const T _xz,
	      const T _yx, const T _yy, const T _yz,
	      const T _zx, const T _zy, const T _zz)
	: xx(_xx), xy(_xy), xz(_xz),
	  yx(_yx), yy(_yy), yz(_yz),
	  zx(_zx), zy(_zy), zz(_zz){}
      Matrix3(const T s) : xx(s), xy(s), xz(s),
			   yx(s), yy(s), yz(s),
			   zx(s), zy(s), zz(s) {}
      Matrix3(const Matrix3 & src)
	: xx(src.xx), xy(src.xy), xz(src.xz),
	  yx(src.yx), yy(src.yy), yz(src.yz),
	  zx(src.zx), zy(src.zy), zz(src.zz) {}

        const Matrix3 & operator = (const Matrix3 & rhs) {
            xx = rhs.xx; xy = rhs.xy; xz = rhs.xz;
	    yx = rhs.yx; yy = rhs.yy; yz = rhs.yz;
	    zx = rhs.zx; zy = rhs.zy; zz = rhs.zz;
            return (*this);
        }
        const Matrix3 & operator = (const T s) {
	  xx = yx = zx = s;
	  xy = yy = zy = s;
	  xz = yz = zz = s;
	  return (*this);
        }

        Matrix3 operator + (const Matrix3 & rhs) const {
	  return Matrix3(xx + rhs.xx, xy + rhs.xy, xz + rhs.xz,
			 yx + rhs.yx, yy + rhs.yy, yz + rhs.yz,
			 zx + rhs.zx, zy + rhs.zy, zz + rhs.zz);
        }
      /*
	Matrix3 operator + (const T & rhs) const {
	  return Matrix3(xx + rhs, xy + rhs, xz + rhs,
			 yx + rhs, yy + rhs, yz + rhs,
			 zx + rhs, zy + rhs, zz + rhs)

        }
      Matrix3 operator - (const T & rhs) const {
	return Matrix3(xx - rhs, yy - rhs, xy, yx);
      }
      const Matrix3 & operator += (const Matrix3 & rhs) {
	(*this) = (*this) + rhs;
	return (*this);
      }
        Matrix3 operator - (const Matrix3 & rhs) const {
            return Matrix3(xx - rhs.xx, yy - rhs.yy, xy - rhs.xy, yx - rhs.yx);
        }
        const Matrix3 & operator -= (const Matrix3 & rhs) {
            (*this) = (*this) - rhs;
            return (*this);
        }
      //*/
      Vector3<T> operator * (const Vector3<T> & rhs) const {
	return Vector3<T>(Vector3<T>(xx,xy,xz)*rhs,
			  Vector3<T>(yx,yy,yz)*rhs,
			  Vector3<T>(zx,zy,zz)*rhs);
      }
      Matrix3 operator * (const Matrix3 & rhs) const {
	return Matrix3(Vector3<T>(xx,xy,xz)*Vector3<T>(rhs.xx,rhs.yx,rhs.zx),
		       Vector3<T>(xx,xy,xz)*Vector3<T>(rhs.xy,rhs.yy,rhs.zy),
		       Vector3<T>(xx,xy,xz)*Vector3<T>(rhs.xz,rhs.yz,rhs.zz),

		       Vector3<T>(yx,yy,yz)*Vector3<T>(rhs.xx,rhs.yx,rhs.zx),
		       Vector3<T>(yx,yy,yz)*Vector3<T>(rhs.xy,rhs.yy,rhs.zy),
		       Vector3<T>(yx,yy,yz)*Vector3<T>(rhs.xz,rhs.yz,rhs.zz),

		       Vector3<T>(zx,zy,zz)*Vector3<T>(rhs.xx,rhs.yx,rhs.zx),
		       Vector3<T>(zx,zy,zz)*Vector3<T>(rhs.xy,rhs.yy,rhs.zy),
		       Vector3<T>(zx,zy,zz)*Vector3<T>(rhs.xz,rhs.yz,rhs.zz));
      }

      const Matrix3 operator - () const {
      return Matrix3(-xx, -xy, -xz,
		     -yx, -yy, -yz,
		     -zx, -zy, -zz);
      }


      friend T Trace(const Matrix3& m){
      return (m.xx + m.yy + m.zz);
      }

      friend const Matrix3 Trans(const Matrix3& m){
	return Matrix3(m.xx, m.yx, m.zx,
		       m.xy, m.yy, m.zy,
		       m.xz, m.yz, m.zz);
      }
    };
  typedef Matrix3<double> F64mat3asym;
}
