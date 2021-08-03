#pragma once

#include"matrix3.hpp"
#include"my_lib.hpp"

namespace MY_LIB{
    // a: semi-major axis
    // l: mean anomaly
    // e: eccentricity
    // u: eccentric anomaly. pericenter for u=0, apocenter for u=pi
    // n: mean mortion
    
    template<typename T>
    T keplereq(const T l,  const T e, const T u){
	return (u - e*sin(u) - l);
    }
    
    template<typename T>
    T keplereq_dot(const T e, const T u){
	return ( 1.0 - e*cos(u) );
    }

    // return eccentric anomaly (u)
    template<typename T>
    T solve_keplereq(const T l,
		     const T e){
	T u0 = l;
	T u1;
	int loop = 0;
	while(1){
	    loop++;
	    T su0 = sin(u0);
	    T cu0 = cos(u0);
	    //T cu0 = sqrt(1.0 - su0*su0);
	    u1 = u0 - ((u0-e*su0-l)/(1.0 - e*cu0));
	    if( fabs(u1-u0) < 1e-15 ){ return u1; }
	    else{ u0 = u1; }
	}
    }

    // Center of mass is located at the origin.
    template<typename Tvec, typename T>
    void OrbParam2PosVel(Tvec & pos0,   Tvec & pos1,
			 Tvec & vel0,   Tvec & vel1,
			 const T mass0, const T mass1,
			 const T ax,    const T ecc,
			 const T inc,   const T OMG,
			 const T omg,   const T u){
	const auto m_tot = mass0 + mass1;
	const auto n = sqrt( m_tot / (ax*ax*ax) );
	const auto cosu = cos(u);
	const auto sinu = sin(u);
	const auto c0 = sqrt(1.0 - ecc*ecc);
	Tvec pos_star(ax*(cosu - ecc), ax*c0*sinu, 0.0);
	Tvec vel_star(-ax*n*sinu/(1.0-ecc*cosu), ax*n*c0*cosu/(1.0-ecc*cosu), 0.0);
	Matrix3<T> rot;
	rot.rotation_zxz(OMG, inc, omg);
	Tvec pos_red = rot*pos_star;
	Tvec vel_red = rot*vel_star;
	pos0 = - mass1 / m_tot * pos_red;
	pos1 =  mass0 / m_tot * pos_red;
	vel0 = - mass1 / m_tot * vel_red;
	vel1 =  mass0 / m_tot * vel_red;
    }

    
}
