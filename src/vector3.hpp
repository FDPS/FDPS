#pragma once

#include<iostream>
#include<iomanip>

namespace ParticleSimulator{
    template<class T>
    class Vector3{
    public:
        T x, y, z;
        Vector3() : x(T(0)), y(T(0)), z(T(0)) {}
        Vector3(const T _x, const T _y, const T _z) : x(_x), y(_y), z(_z) {}
        Vector3(const T s) : x(s), y(s), z(s) {}
        Vector3(const Vector3 & src) : x(src.x), y(src.y), z(src.z) {}

        typedef T DataType;
        static const int DIM = 3;

        const Vector3 & operator = (const Vector3 & rhs){
            x = rhs.x;
            y = rhs.y;
            z = rhs.z;
            return (*this);
        }

        const Vector3 & operator = (const T s){
            x = y = z = s;
            return (*this);
        }

        Vector3 operator + (const Vector3 & rhs) const{
            return Vector3(x + rhs.x, y + rhs.y, z + rhs.z);
        }
        const Vector3 & operator += (const Vector3 & rhs){
            (*this) = (*this) + rhs;
            return (*this);
        }
        Vector3 operator - (const Vector3 & rhs) const{
            return Vector3(x - rhs.x, y - rhs.y, z - rhs.z);
        }

        const Vector3 & operator -= (const Vector3 & rhs){
            (*this) = (*this) - rhs;
            return (*this);
        }

        Vector3 operator * (const T s) const{
            return Vector3(x * s, y * s, z * s);
        }
        const Vector3 & operator *= (const T s){
            (*this) = (*this) * s;
            return (*this);
        }
        friend Vector3 operator * (const T s, const Vector3 & v){
            return (v * s);
        }
        Vector3 operator / (const T s) const{
            return Vector3(x / s, y / s, z / s);
        }
        const Vector3 & operator /= (const T s){
            (*this) = (*this) / s;
            return (*this);
        }

        const Vector3 & operator + () const {
            return (* this);
        }

        const Vector3 operator - () const {
            return Vector3(-x, -y, -z);
        }

        T operator * (const Vector3 & rhs) const{
            return (x * rhs.x) + (y * rhs.y) + (z * rhs.z);
        }

        Vector3 operator ^ (const Vector3 & rhs) const{
            return Vector3( (y * rhs.z - z * rhs.y),
                            (z * rhs.x - x * rhs.z),
                            (x * rhs.y - y * rhs.x) );
        }

        template <typename U>
        operator Vector3<U> () const {
            return Vector3<U> (static_cast<U>(x),
                               static_cast<U>(y),
                               static_cast<U>(z));
        }

        T getMax() const {
            T max_val = (x > y) ? x : y;
            max_val = (max_val > z ) ? max_val : z;
            return max_val;
        }

        T getMin() const {
            T min_val = (x < y) ? x : y;
            min_val = (min_val < z ) ? min_val : z;
            return min_val;
        }

        template <class F>
        Vector3 applyEach(F f) const {
            return Vector3(f(x), f(y), f(z));
        }

        template <class F>
        friend Vector3 ApplyEach(F f, const Vector3 & arg1, const Vector3 & arg2){
            return Vector3( f(arg1.x, arg2.x), f(arg1.y, arg2.y), f(arg1.z, arg2.z) );
        }

        friend std::ostream & operator <<(std::ostream & c, const Vector3 & u){
            c<<u.x<<"   "<<u.y<<"    "<<u.z;
            return c;
        }

	friend std::istream & operator >>(std::istream & c, Vector3 & u){
            c>>u.x; c>>u.y; c>>u.z;
            return c;
        }

#if 0
        const T & operator[](const int i) const {
#ifdef PARTICLE_SIMULATOR_VECTOR_RANGE_CHECK
	    if(i >= DIM || i < 0){
		std::cout<<"PS_ERROR: Vector invalid access. \n"<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;		
		std::cerr<<"Vector element="<<i<<" is not valid."<<std::endl;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
		MPI_Abort(MPI_COMM_WORLD,-1);
#else
		exit(-1);
#endif	//PARTICLE_SIMULATOR_MPI_PARALLEL
	    }
#endif //PARTICLE_SIMULATOR_VECTOR_RANGE_CHECK
            return (&x)[i];
        }

        T & operator[](const int i){
#ifdef PARTICLE_SIMULATOR_VECTOR_RANGE_CHECK
	    if(i >= DIM || i < 0){
		std::cout<<"PS_ERROR: Vector invalid access. \n"<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;		
		std::cerr<<"Vector element="<<i<<" is not valid."<<std::endl;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
		MPI_Abort(MPI_COMM_WORLD,-1);
#else
		exit(-1);
#endif //PARTICLE_SIMULATOR_MPI_PARALLEL		
	    }
#endif //PARTICLE_SIMULATOR_VECTOR_RANGE_CHECK
            return (&x)[i];
        }
#else // #if 0
        const T & operator[](const int i) const {
            if(0==i) return x;
            if(1==i) return y;
            if(2==i) return z;
            std::cout<<"PS_ERROR: Vector invalid access. \n"<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;		
            std::cerr<<"Vector element="<<i<<" is not valid."<<std::endl;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            MPI_Abort(MPI_COMM_WORLD,-1);
#endif //PARTICLE_SIMULATOR_MPI_PARALLEL
            exit(-1);
            return x; //dummy for avoid warning
        }

        T & operator[](const int i){
            //std::cerr<<"operator []"<<std::endl;
            if(0==i) return x;
            if(1==i) return y;
            if(2==i) return z;
            std::cout<<"PS_ERROR: Vector invalid access. \n"<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;		
            std::cerr<<"Vector element="<<i<<" is not valid."<<std::endl;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            MPI_Abort(MPI_COMM_WORLD,-1);
#endif //PARTICLE_SIMULATOR_MPI_PARALLEL
            exit(-1);
            return x; //dummy for avoid warning
        }
#endif // #if 0

        T getDistanceSQ(const Vector3 & u) const {
            T dx = x - u.x;
            T dy = y - u.y;
            T dz = z - u.z;
            return dx*dx + dy*dy + dz*dz;
        }
        bool operator == (const Vector3 & u) const {
            return ( (x==u.x) && (y==u.y) && (z==u.z) );
        }
        bool operator != (const Vector3 & u) const {
            return ( (x!=u.x) || (y!=u.y) || (z!=u.z) );
        }
	
	/*
	Vector3 getDiagonal (const Vector3 & u) const {
	    return Vector3(x*u.x, y*u.y, z*u.z);
	}
	Vector3 getDiagonal (const Vector3<int> & u) const {
	    return Vector3(x*(T)(u.x), y*(T)(u.y), z*(T)(u.z));
	}
	*/
    };

    template <>
    inline Vector3<float> Vector3<float>::operator / (const float s) const {
        const float inv_s = 1.0f/s;
        return Vector3(x * inv_s, y * inv_s, z * inv_s);
    }
    template <>
    inline Vector3<double> Vector3<double>::operator / (const double s) const {
        const double inv_s = 1.0/s;
        return Vector3(x * inv_s, y * inv_s, z * inv_s);
    }
}
