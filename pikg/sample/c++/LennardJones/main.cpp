#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <random>
#include <cmath>
#include <cassert>
#include <chrono>

struct Particle{
  double rx,ry,rz;
  double vx,vy,vz;
  double fx,fy,fz;
  double p;
};

struct EP{
  double rx,ry,rz;
};

struct Force{
  double fx,fy,fz;
  double p;
};

struct Energy{
  double p,k;
};

void clear(Force* f,const int ni){
  for(int i=0;i<ni;i++){
    f[i].fx = 0.0;
    f[i].fy = 0.0;
    f[i].fz = 0.0;
    f[i].p  = 0.0;
  }
}

void copyToEp(EP* ep,const Particle* ptcls,const int ni){
  for(int i=0;i<ni;i++){
    ep[i].rx = ptcls[i].rx;
    ep[i].ry = ptcls[i].ry;
    ep[i].rz = ptcls[i].rz;
  }
}

void copyFromForce(Particle* ptcls,const Force* f,const int ni){
  for(int i=0;i<ni;i++){
    ptcls[i].fx = f[i].fx;
    ptcls[i].fy = f[i].fy;
    ptcls[i].fz = f[i].fz;
    ptcls[i].p  = f[i].p;
  }
}

Energy calcEnergy(const Particle* ptcls,const int N){
  Energy e;
  e.p = e.k = 0.0;
  for(int i=0;i<N;i++){
    e.p += ptcls[i].p;
    e.k += ptcls[i].vx * ptcls[i].vx
         + ptcls[i].vy * ptcls[i].vy
         + ptcls[i].vz * ptcls[i].vz;
  }
  e.p *= 0.5;
  e.k *= 0.5;
  return e;
}

void scaleVelocity(Particle* ptcls,const int N,const double T){
  const double ktmp = calcEnergy(ptcls,N).k;
  const double scale = sqrt(1.5*N*T/ktmp);
  double vtotx,vtoty,vtotz;
  vtotx = vtoty = vtotz = 0.0;
  for(int i=0;i<N;i++){
    ptcls[i].vx *= scale;
    ptcls[i].vy *= scale;
    ptcls[i].vz *= scale;
    vtotx += ptcls[i].vx;
    vtoty += ptcls[i].vy;
    vtotz += ptcls[i].vz;
  }
  // cancel total momentum
  vtotx /= (double)N;
  vtoty /= (double)N;
  vtotz /= (double)N;
  for(int i=0;i<N;i++){
    ptcls[i].vx -= vtotx;
    ptcls[i].vy -= vtoty;
    ptcls[i].vz -= vtotz;
  }
}

#include "kernel.hpp"

int main(int argc,char **argv){
  const int n = 10;
  const int N = n*n*n;
  const double rho = 1.05;
  const double l = std::pow((double)N/rho,1./3.);
  const double lh = 0.5*l;
  const double u = l / (double)n;
  const double T = 2.0;
  const double rc = lh > 4.5 ? 4.5 : lh;

  std::cerr << "N= " << N << std::endl;
  std::cerr << "l= " << l << " or " << "rho= " << rho << std::endl;
  std::cerr << "T= " << T << std::endl;

  Particle* ptcls = new Particle[N];
  EP* ep = new EP[N];
  Force* force = new Force[N];

  int count = 0;
  std::random_device seed_gen;
  std::mt19937 engine(seed_gen());
  std::uniform_real_distribution<> dist(-1.0, 1.0);

  for(int x=0;x<n;x++){
    for(int y=0;y<n;y++){
      for(int z=0;z<n;z++){
	ptcls[count].rx = u*x - lh;
	ptcls[count].ry = u*y - lh;
	ptcls[count].rz = u*z - lh;

	ptcls[count].vx = dist(engine);
	ptcls[count].vy = dist(engine);
	ptcls[count].vz = dist(engine);

	ptcls[count].fx = 0.0;
	ptcls[count].fy = 0.0;
	ptcls[count].fz = 0.0;
	ptcls[count].p  = 0.0;
	count++;
      }
    }
  }
  assert(count == N);
  scaleVelocity(ptcls,N,T);

  copyToEp(ep,ptcls,N);
  clear(force,N);
  Kernel calcForce(l,rc);
  calcForce(ep,N,ep,N,force);
  copyFromForce(ptcls,force,N);

  const double dt = 0.0005;
  const double dth = 0.5*dt;
  std::cerr << "dt= " << dt << std::endl;
  std::cerr << "dth= " << dth << std::endl;
  const int nstep = 2000;
  double etot = 0.0;
  double time = 0.0;
  for(int s=-2000;s<nstep;s++){
    //std::cerr << s << std::endl;
    // equilibration
    if(s < 0) scaleVelocity(ptcls,N,T);
    // first kick
    for(int i=0;i<N;i++){
      ptcls[i].vx += ptcls[i].fx * dth;
      ptcls[i].vy += ptcls[i].fy * dth;
      ptcls[i].vz += ptcls[i].fz * dth;
      //std::cout << ptcls[i].fx << " " << ptcls[i].fy << " " << ptcls[i].fz << std::endl;
    }
    // drift
    //std::cout << "pbc" << std::endl;
    for(int i=0;i<N;i++){
      ptcls[i].rx += ptcls[i].vx * dt;
      ptcls[i].ry += ptcls[i].vy * dt;
      ptcls[i].rz += ptcls[i].vz * dt;
      //std::cout << ptcls[i].rx << " " << ptcls[i].ry << " " << ptcls[i].rz << std::endl;
      // periodic boundary condition
      while(ptcls[i].rx < -lh) ptcls[i].rx += l;
      while(ptcls[i].rx >= lh) ptcls[i].rx -= l;
      while(ptcls[i].ry < -lh) ptcls[i].ry += l;
      while(ptcls[i].ry >= lh) ptcls[i].ry -= l;
      while(ptcls[i].rz < -lh) ptcls[i].rz += l;
      while(ptcls[i].rz >= lh) ptcls[i].rz -= l;
    }
    //std::cout << "calcForce" << std::endl;
    // calc force and pot
    copyToEp(ep,ptcls,N);
    clear(force,N);
    std::chrono::system_clock::time_point start,end;
    start = std::chrono::system_clock::now();
    calcForce(ep,N,ep,N,force);
    end = std::chrono::system_clock::now();
    if(s>=0) time += std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
    copyFromForce(ptcls,force,N);

    // second kick
    for(int i=0;i<N;i++){
      ptcls[i].vx += ptcls[i].fx * dth;
      ptcls[i].vy += ptcls[i].fy * dth;
      ptcls[i].vz += ptcls[i].fz * dth;
    }

    if(s%100 == 0 && s>=0){
      Energy e = calcEnergy(ptcls,N);
      const double etmp = e.p + e.k;
      if(s==0) etot = etmp;
      std::cout << std::scientific << s*dt << " " << e.p << " " << e.k << " " << etmp << " " << std::abs((etmp - etot)/etot) << std::endl;
#if 0
      // output cdv
      std::stringstream strs;
      strs << "test" << std::setw(8) << std::setfill('0') << s/100 << ".cdv";
      std::ofstream ofs(strs.str());
      ofs << "'box_sx=" << -lh << ",box_ex=" << lh << std::endl;
      ofs << "'box_sy=" << -lh << ",box_ey=" << lh << std::endl;
      ofs << "'box_sz=" << -lh << ",box_ez=" << lh << std::endl;
      for(int i=0;i<N;i++) ofs << i << " 0 " << ptcls[i].rx << " " << ptcls[i].ry << " " << ptcls[i].rz << std::endl;
#endif
    }
  }
  // check total energy error
  Energy e = calcEnergy(ptcls,N);
  const double error = std::abs((e.p + e.k - etot)/etot);
  if(error > 1.5e-4){ std::cerr << "test failed" << std::endl; return -1;}
  else                std::cerr << "test passed" << std::endl;

  // print flops
  //std::cerr << "elapsed_time_for_one_kernel(ns): " << time/nstep << std::endl;
  std::cerr << "kernel_performance(GFlops): " << 26.0 * N * N * nstep / (time*1e-9) / 1e9 << std::endl;
  return 0;
}
