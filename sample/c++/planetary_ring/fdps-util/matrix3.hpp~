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
	
	/*
	template<typename Tvec>
	friend Matrix3 outer_product(const Tvec & a, const Tvec &b){
	    Matrix3 x;
	    x[0][0]=a[0]*b[0]; x[0][1]=a[0]*b[1]; x[0][2]=a[0]*b[2];
	    x[1][0]=a[1]*b[0]; x[1][1]=a[1]*b[1]; x[1][2]=a[1]*b[2];
	    x[2][0]=a[2]*b[0]; x[2][1]=a[2]*b[1]; x[2][2]=a[2]*b[2];
	    return x;
	}
	*/
    /*
    void rotation(const T I, const T OMEGA, const T omega){
	const T cosomega = cos(omega);
	const T sinomega = sin(omega);
	const T cosOMEGA = cos(OMEGA);
	const T sinOMEGA = sin(OMEGA);
	const T cosinc = cos(I);
	const T sininc = sin(I);
	xx = cosomega*cosOMEGA  - sinomega*sinOMEGA*cosinc;
	xy = -sinomega*cosOMEGA - cosomega*sinOMEGA*cosinc;
	xz = sinOMEGA*sininc;

	yx = cosomega*sinOMEGA  + sinomega*cosOMEGA*cosinc;
	yy = -sinomega*sinOMEGA + cosomega*cosOMEGA*cosinc;
	yz = -cosOMEGA*sininc;

	zx = sinomega*sininc;
	zy = cosomega*sininc;
	zz = cosinc;
    }
    */


    /*

    */



    /*
    const T * operator[](const int & i)const{return m[i];}
    
    T * operator[](const int & i){return m[i];}

    Matrix3& operator=(const Matrix3 &c)const{
	m[0][0] = c.m[0][0]; m[0][1] = c.m[0][1]; m[0][2] = c.m[0][2];
	m[1][0] = c.m[1][0]; m[1][1] = c.m[1][1]; m[1][2] = c.m[1][2];
	m[2][0] = c.m[2][0]; m[2][1] = c.m[2][1]; m[2][2] = c.m[2][2];
	return *this;
    }

    Matrix3 & operator+(const Matrix3 &a)const{
	Matrix3 x;
	x.m[0][0]=m[0][0]+a.m[0][0]; x.m[0][1]=m[0][1]+a.m[0][1];  x.m[0][2]=m[0][2]+a.m[0][2];
	x.m[1][0]=m[1][0]+a.m[1][0]; x.m[1][1]=m[1][1]+a.m[1][1];  x.m[1][2]=m[1][2]+a.m[1][2];
	x.m[2][0]=m[2][0]+a.m[2][0]; x.m[2][1]=m[2][1]+a.m[2][1];  x.m[2][2]=m[2][2]+a.m[2][2];
	return x;
    }

    matrix3& operator-(const matrix3 &a)const{
	matrix3 x;
	x.m[0][0]=m[0][0]-a.m[0][0]; x.m[0][1]=m[0][1]-a.m[0][1];  x.m[0][2]=m[0][2]-a.m[0][2];
	x.m[1][0]=m[1][0]-a.m[1][0]; x.m[1][1]=m[1][1]-a.m[1][1];  x.m[1][2]=m[1][2]-a.m[1][2];
	x.m[2][0]=m[2][0]-a.m[2][0]; x.m[2][1]=m[2][1]-a.m[2][1];  x.m[2][2]=m[2][2]-a.m[2][2];
	return x;
    }


    matrix3& operator*(const matrix3 &a)const{
	matrix3 x;
	x.m[0][0] = m[0][0]*a.m[0][0] + m[0][1]*a.m[1][0] + m[0][2]*a.m[2][0];
	x.m[0][1] = m[0][0]*a.m[0][1] + m[0][1]*a.m[1][1] + m[0][2]*a.m[2][1];
	x.m[0][2] = m[0][0]*a.m[0][2] + m[0][1]*a.m[1][2] + m[0][2]*a.m[2][2];

	x.m[1][0] = m[1][0]*a.m[0][0] + m[1][1]*a.m[1][0] + m[1][2]*a.m[2][0];
	x.m[1][1] = m[1][0]*a.m[0][1] + m[1][1]*a.m[1][1] + m[1][2]*a.m[2][1];
	x.m[1][2] = m[1][0]*a.m[0][2] + m[1][1]*a.m[1][2] + m[1][2]*a.m[2][2];

	x.m[2][0] = m[2][0]*a.m[0][0] + m[2][1]*a.m[1][0] + m[2][2]*a.m[2][0];
	x.m[2][1] = m[2][0]*a.m[0][1] + m[2][1]*a.m[1][1] + m[2][2]*a.m[2][1];
	x.m[2][2] = m[2][0]*a.m[0][2] + m[2][1]*a.m[1][2] + m[2][2]*a.m[2][2];
	return x;
    }


    vector3<T>& operator * (const vector3<T> &a)const{
	vector3<T> x;
	x[0] = m[0][0]*a[0] + m[0][1]*a[1] + m[0][2]*a[2];
	x[1] = m[1][0]*a[0] + m[1][1]*a[1] + m[1][2]*a[2];
	x[2] = m[2][0]*a[0] + m[2][1]*a[1] + m[2][2]*a[2];
	return x;
    }

    matrix3& operator * (const T &s) const{
	matrix3 x;
	x.m[0][0]=m[0][0]*s; x.m[0][1]=m[0][1]*s; x.m[0][2]=m[0][2]*s;
	x.m[1][0]=m[1][0]*s; x.m[1][1]=m[1][1]*s; x.m[1][2]=m[1][2]*s;
	x.m[2][0]=m[2][0]*s; x.m[2][1]=m[2][1]*s; x.m[2][2]=m[2][2]*s;
	return x;
    }

    matrix3& transposed()const{
	matrix3 mtmp;
	for(int i=0; i<3; i++){
	    for(int j=0; j<3; j++){
		mtmp.m[j][i] = m[i][j];
	    }
	}
	return mtmp;
    }

    void unit(){
	m[0][1] = m[0][2] = m[1][0] = m[1][2] = m[2][0] = m[2][1] = 0.0;
	m[0][0] = m[1][1] = m[2][2] = 1.0;
    }

    matrix3& operator += (const matrix3 &a){
	m[0][0] += a.m[0][0]; m[0][1] += a.m[0][1];  m[0][2] += a.m[0][2];
	m[1][0] += a.m[1][0]; m[1][1] += a.m[1][1];  m[1][2] += a.m[1][2];
	m[2][0] += a.m[2][0]; m[2][1] += a.m[2][1];  m[2][2] +=a.m[2][2];
	return *this;
    }

    matrix3& operator -= (const matrix3 &a){
	m[0][0] -= a.m[0][0]; m[0][1] -= a.m[0][1];  m[0][2] -= a.m[0][2];
	m[1][0] -= a.m[1][0]; m[1][1] -= a.m[1][1];  m[1][2] -= a.m[1][2];
	m[2][0] -= a.m[2][0]; m[2][1] -= a.m[2][1];  m[2][2] -=a.m[2][2];
	return *this;
    }

    matrix3& operator *= (const T &a){
	m[0][0] *= a; m[0][1] *= a;  m[0][2] *= a;
	m[1][0] *= a; m[1][1] *= a;  m[1][2] *= a;
	m[2][0] *= a; m[2][1] *= a;  m[2][2] *= a;
	return *this;
    }

    T& determinant()const{
	T det = 0.0;
	        det = m[0][0] * m[1][1] * m[2][2]
		                + m[0][1] * m[1][2] * m[2][0]
		                + m[0][2] * m[1][0] * m[2][1]
		                - m[0][2] * m[1][1] * m[2][0]
		                - m[0][1] * m[1][0] * m[2][2]
		    - m[0][0] * m[1][2] * m[2][1];
		return det;
    }



    friend matrix3 quadrupole(const vector3<T> &a){
	matrix3 x=3.0*outer_product(a,a);
	T a_squared=a*a;
	matrix3 delta;
	delta.unit();
	return x-a_squared*delta;
    }

    friend matrix3 inertiamoment(const T mass, const vector3<T> &pos){
	matrix3 x = mass * outer_product(pos, pos);
	return x;
    }

    friend vector3<T> operator*(const vector3<T> &a, const matrix3 &c){
	vector3<T> x;
	x[0]=a[0]*c.m[0][0] + a[1]*c.m[1][0] + a[2]*c.m[2][0];
	x[1]=a[0]*c.m[0][1] + a[1]*c.m[1][1] + a[2]*c.m[2][1];
	x[2]=a[0]*c.m[0][2] + a[1]*c.m[1][2] + a[2]*c.m[2][2];
	return x;
    }

    friend matrix3 operator * (const T &s, const matrix3 &c){
	matrix3 x;
	x[0][0]=c.m[0][0]*s; x[0][1]=c.m[0][1]*s; x[0][2]=c.m[0][2]*s;
	x[1][0]=c.m[1][0]*s; x[1][1]=c.m[1][1]*s; x[1][2]=c.m[1][2]*s;
	x[2][0]=c.m[2][0]*s; x[2][1]=c.m[2][1]*s; x[2][2]=c.m[2][2]*s;
	return x;
    }
    */    
    };
    
}
