#pragma once

#include<iostream>
#include<iomanip>

namespace ParticleSimulator{

    template <typename T>
    class Vector2{
    public:
        T x, y;
        Vector2() : x(T(0)), y(T(0)) {}
        Vector2(const T _x, const T _y) : x(_x), y(_y) {}
        Vector2(const T s) : x(s), y(s) {}
        Vector2(const Vector2 & src) : x(src.x), y(src.y) {}

        static const int DIM = 2;
	
        const Vector2 & operator = (const Vector2 & rhs){
            x = rhs.x;
            y = rhs.y;
            return (*this);
        }

        const Vector2 & operator = (const T s){
            x = y = s;
            return (*this);
        }

        Vector2 operator + (const Vector2 & rhs) const{
            return Vector2(x + rhs.x, y + rhs.y);
        }
        const Vector2 & operator += (const Vector2 & rhs){
            (*this) = (*this) + rhs;
            return (*this);
        }
        Vector2 operator - (const Vector2 & rhs) const{
            return Vector2(x - rhs.x, y - rhs.y);
        }
        const Vector2 & operator -= (const Vector2 & rhs){
            (*this) = (*this) - rhs;
            return (*this);
        }

        // vector scholar products
        Vector2 operator * (const T s) const{
            return Vector2(x * s, y * s);
        }
        const Vector2 & operator *= (const T s){
            (*this) = (*this) * s;
            return (*this);
        }
        friend Vector2 operator * (const T s, const Vector2 & v){
            return (v * s);
        }
        Vector2 operator / (const T s) const{
            return Vector2(x / s, y / s);
        }
        const Vector2 & operator /= (const T s){
            (*this) = (*this) / s;
            return (*this);
        }

        const Vector2 & operator + () const {
            return (* this);
        }

        const Vector2 operator - () const {
            return Vector2(-x, -y);
        }

        // inner product
        T operator * (const Vector2 & rhs) const{
            return (x * rhs.x) + (y * rhs.y);
        }

        // outer product (retruned value is scholar)
        T operator ^ (const Vector2 & rhs) const{
            const T z = (x * rhs.y) - (y * rhs.x);
            return z;
        }

        //cast to Vector2<U>
        template <typename U>
        operator Vector2<U> () const {
            return Vector2<U> (static_cast<U>(x),
                               static_cast<U>(y));
        }

        T getMax() const {
            return x > y ? x : y;
        }

        T getMin() const {
            return x < y ? x : y;
        }

        template <class F>
        Vector2 applyEach(F f) const {
            return Vector2(f(x), f(y));
        }
        template <class F>
        friend Vector2 ApplyEach(F f, const Vector2 & arg1, const Vector2 & arg2){
            return Vector2( f(arg1.x, arg2.x), f(arg1.y, arg2.y) );
        }

        friend std::ostream & operator <<(std::ostream & c, const Vector2 & u){
            c<<u.x<<"   "<<u.y;
            return c;
        }

	friend std::istream & operator >>(std::istream & c, Vector2 & u){
            c>>u.x; c>>u.y;
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
#endif		
	    }
#endif
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
#endif		
	    }	    
#endif	    
            return (&x)[i];
        }
#else
        const T & operator[](const int i) const {
			if(0==i) return x;
			if(1==i) return y;
			std::cout<<"PS_ERROR: Vector invalid access. \n"<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;		
			std::cerr<<"Vector element="<<i<<" is not valid."<<std::endl;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
			MPI_Abort(MPI_COMM_WORLD,-1);
#endif		
			exit(-1);
            return x; //dummy for avoid warning
		}
        T & operator[](const int i){
			if(0==i) return x;
			if(1==i) return y;
			std::cout<<"PS_ERROR: Vector invalid access. \n"<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;		
			std::cerr<<"Vector element="<<i<<" is not valid."<<std::endl;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
			MPI_Abort(MPI_COMM_WORLD,-1);
#endif		
			exit(-1);
            return x; //dummy for avoid warning
		}
#endif

	T getDistanceSQ(const Vector2 & u) const {
            T dx = x - u.x;
            T dy = y - u.y;
            return dx*dx + dy*dy;
        }
        bool operator == (const Vector2 & u) const {
            return ( (x==u.x) && (y==u.y) );
        }
        bool operator != (const Vector2 & u) const {
            return ( (x!=u.x) || (y!=u.y) );
        }

    };

    template <>
    inline Vector2<float> Vector2<float>::operator / (const float s) const {
        const float inv_s = 1.0f/s;
        return Vector2(x * inv_s, y * inv_s);
    }
    template <>
    inline Vector2<double> Vector2<double>::operator / (const double s) const {
        const double inv_s = 1.0/s;
        return Vector2(x * inv_s, y * inv_s);
    }
}
