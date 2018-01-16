#ifndef H_UNIT
#define H_UNIT
const double Na  = 6.0221367e+23;     // avogadro number                / -
const double e0  = 8.8542e-12;        // permeability of vacuum         / A^2 s^4 / kg m^3
const double ele = 1.6021761e-19;     // elementary charge (denkisoryo) / C
const double k_B = 1.3806488e-23;     // boltzmann constant             / J / K
const double k_e = 1.0/(4.0*M_PI*e0); // coulomb's constant             / N m^2 / C^2

//const double unit_mass   = 1.00794e-3/Na;   // kg      (mass of hydrogen)
const double unit_mass   = 15.9994e-3/Na;   // kg      (mass of oxygen)
const double unit_length = 1.0e-10;         // m       (1 angstrom)
const double unit_energy = 1e3 / Na;        // J / mol (1 kJ/mol)

const double unit_density= unit_mass/(unit_length*unit_length*unit_length);
const double unit_time   = sqrt(unit_mass*unit_length*unit_length / unit_energy); // s
const double unit_coulomb= sqrt(4.*M_PI*e0*unit_energy*unit_length);  // sqrt(C^2 kJ/mol)
const double unit_temp   = unit_energy / k_B;                                     // K
const double unit_press  = unit_energy / (unit_length*unit_length*unit_length);   // Pa
#endif
