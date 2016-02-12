#ifndef __GADGET_PARAM__
#define __GADGET_PARAM__

namespace ParticleSimulator{
    namespace ParticleMesh{

//#define LONGIDS

const double Gadget_UnitLength_in_Mpc = 0.001;          //  1kpc
//const double Gadget_UnitLength_in_Mpc = 1.0;            // 1Mpc
const double Gadget_UnitMass_in_Msun = 1.0e10;          // 1e10 Msun
const double Gadget_UnitVelocity_in_cm_per_s = 1e5;     //  1 km/sec


typedef struct GadgetHeader{

  int      Npart[6];
  double   Massarr[6];
  double   Time;
  double   Redshift;
  int      FlagSfr;
  int      FlagFeedback;
  unsigned int Nall[6];
  int      FlagCooling;
  int      NumFiles;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam;
  int      Flag_StellarAge;
  int      Flag_Metals;
  unsigned int NallHW[6];
  int      flag_entr_ics;
  char     unused[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8 - 9*4];
  
}GadgetHeader, *pGadgetHeader;



    } // namespace ParticleMesh
}     // namespace ParticleSimulator

#endif
