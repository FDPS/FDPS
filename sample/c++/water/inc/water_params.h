#ifndef H_WATER_PARAMS
#define H_WATER_PARAMS

#include "unit.h"

//SPC/E parameters
const PS::F64 MASS_OXY    = 15.9994 * 1e-3 / Na / unit_mass;
const PS::F64 SIGMA_OXY   = 3.166 * 1e-10 / unit_length;
const PS::F64 EPSILON_OXY = 0.6500 * 1e+3 / Na / unit_energy;
const PS::F64 CHARGE_OXY  = -0.8476 * ele / unit_coulomb;

const PS::F64 MASS_HYD    = 1.00800 * 1e-3 / Na / unit_mass;
const PS::F64 SIGMA_HYD   = 0.0;
const PS::F64 EPSILON_HYD = 0.0;
const PS::F64 CHARGE_HYD  = 0.4238 * ele / unit_coulomb;

const PS::F64 ANGLE_HOH   = M_PI / 180.0 * 109.47;

const PS::F64 BOND_OH = 1.0 * 1e-10 / unit_length;
const PS::F64 BOND_HH = 2.0 * BOND_OH * cos(M_PI/180.0*(90.0 - 0.5*109.47));

const PS::F64 BOND_COEF  = 3450.;
const PS::F64 ANGLE_COEF = 383.;
#endif
