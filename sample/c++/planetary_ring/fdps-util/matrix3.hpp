#pragma once

#include<cmath>
#include<ostream>
#include<iomanip>

namespace MY_LIB{

    template<typename T>
    class Matrix3{
    public:
	T xx, xy, xz, yx, yy, yz, zx, zy, zz;

	void rotation_zxz(const T th0, const T th1, const T th2){
	    const T cos_th2 = cos(th2);
	    const T sin_th2 = sin(th2);
	    const T cos_th0 = cos(th0);
	    const T sin_th0 = sin(th0);
	    const T cos_th1 = cos(th1);
	    const T sin_th1 = sin(th1);

	    xx =  cos_th2*cos_th0 - sin_th2*sin_th0*cos_th1;
	    xy = -sin_th2*cos_th0 - cos_th2*sin_th0*cos_th1;
	    xz =  sin_th0*sin_th1;

	    yx =  cos_th2*sin_th0 + sin_th2*cos_th0*cos_th1;
	    yy = -sin_th2*sin_th0 + cos_th2*cos_th0*cos_th1;
	    yz = -cos_th0*sin_th1;

	    zx = sin_th2*sin_th1;
	    zy = cos_th2*sin_th1;
	    zz = cos_th1;

	    
	}	
	Matrix3() : xx(T(1)), xy(T(0)), xz(T(0)), yx(T(0)), yy(T(1)), yz(T(0)), zx(T(0)), zy(T(0)), zz(T(1)) {}
	Matrix3(const T c) : xx(c), xy(c), xz(c), yx(c), yy(c), yz(c), zx(c), zy(c), zz(c) {}
	Matrix3(const Matrix3 & c) : xx(c.xx), xy(c.xy), xz(c.xz), yx(c.yx), yy(c.yy), yz(c.yz), zx(c.zx), zy(c.zy), zz(c.zz) {}
	Matrix3(const T th0, const T th1, const T th2){
	    rotation_zxz(th0, th1, th2);
	}
	Matrix3(const T _xx, const T _xy, const T _xz,
		const T _yx, const T _yy, const T _yz,
		const T _zx, const T _zy, const T _zz)
	    : xx(_xx), xy(_xy), xz(_xz), yx(_yx), yy(_yy), yz(_yz), zx(_zx), zy(_zy), zz(_zz) {}
	template <typename U>
	operator Matrix3<U> () const {
	    return Matrix3<U>( static_cast<U>(xx), static_cast<U>(xy), static_cast<U>(xz),
			       static_cast<U>(yx), static_cast<U>(yy), static_cast<U>(yz),
			       static_cast<U>(zx), static_cast<U>(zy), static_cast<U>(zz) );
	}

	void euler_rotatoin(const T OMEGA, const T inc, const T omega){
	    rotation_zxz(OMEGA, inc, omega);
	}

	// rotate position itself (not rotate coordinate)
	// Mat * Vec
	template<class Tvec>
	Tvec operator * (const Tvec & vec) const {
	    Tvec ret;
	    ret.x = xx*vec.x + xy*vec.y + xz*vec.z;
	    ret.y = yx*vec.x + yy*vec.y + yz*vec.z;
	    ret.z = zx*vec.x + zy*vec.y + zz*vec.z;
	    return ret;
	}

	friend std::ostream & operator << (std::ostream & c, Matrix3 & mtmp){
	    c<<std::setprecision(15)<<mtmp.xx<<"   "<<mtmp.xy<<"    "<<mtmp.xz<<std::endl;
	    c<<std::setprecision(15)<<mtmp.yx<<"   "<<mtmp.yy<<"    "<<mtmp.yz<<std::endl;
	    c<<std::setprecision(15)<<mtmp.zx<<"   "<<mtmp.zy<<"    "<<mtmp.zz<<std::endl;
	    return c;
	}

	// Vec^T * Mat
	template<class Tvec>
	friend Tvec operator * (const Tvec & vec, const Matrix3 & mat){
	    Tvec ret;
	    ret.x = mat.xx*vec.x + mat.yx*vec.y + mat.zx*vec.z;
	    ret.y = mat.xy*vec.x + mat.yy*vec.y + mat.zy*vec.z;
	    ret.z = mat.xz*vec.x + mat.yz*vec.y + mat.zz*vec.z;
	    return ret;
	}

	void transposed() {
	    std::swap(xy, yx);
	    std::swap(xz, zx);
	    std::swap(yz, zy);
	}
    };
    
}
