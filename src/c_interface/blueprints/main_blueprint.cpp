/* Standard headers */
#include <iostream>
#include <fstream>
/* FDPS headers */
#include <particle_simulator.hpp> 
/* User-defined headers */
#include "FDPS_Manipulators.h"
extern "C" {
int c_main(void);
}
int main(int argc, char *argv[])
{
   
   //* Initialize fdps_manip
   FDPS_Manipulators::Initialize(argc,argv);
   //* Call C main function
   return c_main();

}
