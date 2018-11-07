#include "physical_constants.h"

/* Global variables */
// Physical constants
double Ggrav;  // Gravitational constants [cm^{3}/g/s^{2}]
double kBoltz; // Boltzmann constant [cm^{2} g/s^{2}/K]

// Mass units
double Msolar; // Solar mass [g]
double dalton; // Unified atomic mass unit [g]
double Mhydrogen; // Mass of hydrogen atom [g]

// Length units
double km; // kilo-meters [cm]
double AU; // Astronomical unit [cm]
double pc; // parsec [cm]
double kpc; // kilo-parsecs [cm]
double Mpc; // mega-parsecs [cm]
double Gpc; // giga-parsecs [cm]

// Time units
double yr;   // year [s]
double kyr; // kilo-years [s]
double Myr; // mega-years [s]
double Gyr; // giga-years [s]

void setup_phys_const() {

    // Physical constants
    Ggrav = 6.67408e-8;
    kBoltz = 1.3806503e-16;

    // Mass units
    Msolar = 1.989e33;
    dalton = 1.660538921e-24;
    Mhydrogen = 1.007825e0 * dalton;

    // Length units
    km = 1.0e5;
    AU = 1.49597870e13;
    pc = 3.0857e18;
    kpc = 1.0e3 * pc;
    Mpc = 1.0e6 * pc;
    Gpc = 1.0e9 * pc;

    // Time units
    yr = 3.15576e7;
    kyr = 1.0e3 * yr;
    Myr = 1.0e6 * yr;
    Gyr = 1.0e9 * yr;

}
