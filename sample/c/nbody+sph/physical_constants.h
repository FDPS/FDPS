#pragma once

/* Global variables */
// Physical constants
extern double Ggrav; // Gravitational constants [cm^{3}/g/s^{2}]
extern double kBoltz; // Boltzmann constant [cm^{2} g/s^{2}/K]

// Mass units
extern double Msolar; // Solar mass [g]
extern double dalton; // Unified atomic mass unit [g]
extern double Mhydrogen; // Mass of hydrogen atom [g]

// Length units
extern double km; // kilo-meters [cm]
extern double AU; // Astronomical unit [cm]
extern double pc; // parsec [cm]
extern double kpc; // kilo-parsecs [cm]
extern double Mpc; // mega-parsecs [cm]
extern double Gpc; // giga-parsecs [cm]

// Time units
extern double yr; // year [s]
extern double kyr; // kilo-years [s]
extern double Myr; // mega-years [s]
extern double Gyr; // giga-years [s]


/* Prototype declarations */
void setup_phys_const();
