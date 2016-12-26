extern "C" void f_main_();

/* Standard headers */
#include <iostream>
#include <fstream>
/* FDPS headers */
#include <particle_simulator.hpp> 
/* User-defined headers */
#include "FDPS_Manipulators.h"

int main(int argc, char *argv[])
{
   
   //* Initialize fdps_manip
   FDPS_Manipulators::Initialize(argc,argv);
   //* Call Fortran main subroutine
   f_main_();

   return 0;

}
