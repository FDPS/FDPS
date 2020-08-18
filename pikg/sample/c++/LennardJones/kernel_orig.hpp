struct Kernel{
  const double l;
  const double rc;
  Kernel(const double _l,const double _rc):l(_l),rc(_rc){
    std::cerr << "l= " << l << std::endl;
    std::cerr << "rc= " << rc << std::endl;
    assert(l >= 2.0*rc);
  }
  void operator()(const EP* epi,
		  const int nepi,
		  const EP* epj,
		  const int nepj,
		  Force*    force)
  {
   const double lh = 0.5*l;
   const double rc = 4.0;
   const double rc2 = rc*rc;
   for(int i=0;i<nepi;i++){
     const double xi = epi[i].rx;
     const double yi = epi[i].ry;
     const double zi = epi[i].rz;
     double fx,fy,fz,p;
     fx = fy = fz = p = 0.0;
     for(int j=0;j<nepj;j++){
       double dx = xi - epj[j].rx;
       double dy = yi - epj[j].ry;
       double dz = zi - epj[j].rz;
       if(dx < -lh) dx += l;
       if(dx >= lh) dx -= l;
       if(dy < -lh) dy += l;
       if(dy >= lh) dy -= l;
       if(dz < -lh) dz += l;
       if(dz >= lh) dz -= l;
       const double r2 = dx*dx + dy*dy + dz*dz;

       if(r2 > 0.0 && r2 < rc2){
	 const double r2i = 1.0 / r2;
	 const double r6i = r2i * r2i * r2i;
	 const double f = (48.0*r6i - 24.0)*r6i*r2i;
	 fx += f*dx;
	 fy += f*dy;
	 fz += f*dz;
	 p  += 4.0*r6i*(r6i - 1.0);
       }
     }
     force[i].fx += fx;
     force[i].fy += fy;
     force[i].fz += fz;
     force[i].p  += p;
   }
  }
};
