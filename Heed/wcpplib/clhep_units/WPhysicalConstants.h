// HEP coherent Physical Constants
//
// This file has been provided by Geant4 (simulation toolkit for HEP).
//
// The basic units are :
// millimeter
// nanosecond
// Mega electron Volt
// positon charge
// degree Kelvin
// amount of substance (mole)
// luminous intensity (candela)
// radian
// steradian
//
// Below is a non exhaustive list of Physical CONSTANTS,
// computed in the Internal HEP System Of Units.
//
// Most of them are extracted from the Particle Data Book :
//        Phys. Rev. D  volume 50 3-1 (1994) page 1233
//
//        ...with a meaningful (?) name ...
//
// You can add your own constants.
//
// Author: M.Maire
//
// History:
//
// 23.02.96 Created
// 26.03.96 Added constants for standard conditions of temperature
//          and pressure; also added Gas threshold.

#ifndef HEED_PHYSICAL_CONSTANTS_H
#define HEED_PHYSICAL_CONSTANTS_H

#include "wcpplib/clhep_units/WSystemOfUnits.h"

namespace Heed {

namespace CLHEP {

//
//
//
static constexpr double Avogadro = 6.0221367e+23 / mole;

//
// c   = 299.792458 mm/ns
// c^2 = 898.7404 (mm/ns)^2
//
static constexpr double c_light = 2.99792458e+8 * m / s;
static constexpr double c_squared = c_light * c_light;

//
// h     = 4.13566e-12 MeV*ns
// hbar  = 6.58212e-13 MeV*ns
// hbarc = 197.32705e-12 MeV*mm
//
static constexpr double h_Planck = 6.6260755e-34 * joule * s;
static constexpr double hbar_Planck = h_Planck / twopi;
static constexpr double hbarc = hbar_Planck * c_light;
static constexpr double hbarc_squared = hbarc * hbarc;

//
//
//
static constexpr double electron_charge = -eplus;  // see SystemOfUnits.h
static constexpr double e_squared = eplus * eplus;

//
// amu_c2 - atomic equivalent mass unit
// amu    - atomic mass unit
//
static constexpr double electron_mass_c2 = 0.51099906 * MeV;
static constexpr double proton_mass_c2 = 938.27231 * MeV;
static constexpr double neutron_mass_c2 = 939.56563 * MeV;
static constexpr double amu_c2 = 931.49432 * MeV;
static constexpr double amu = amu_c2 / c_squared;

//
// permeability of free space mu0    = 2.01334e-16 Mev*(ns*eplus)^2/mm
// permittivity of free space epsil0 = 5.52636e+10 eplus^2/(MeV*mm)
//
static constexpr double mu0 = 4 * pi * 1.e-7 * henry / m;
static constexpr double epsilon0 = 1. / (c_squared * mu0);

//
// electromagnetic coupling = 1.43996e-12 MeV*mm/(eplus^2)
//
static constexpr double elm_coupling = e_squared / (4 * pi * epsilon0);
static constexpr double fine_structure_const = elm_coupling / hbarc;
static constexpr double classic_electr_radius = elm_coupling / electron_mass_c2;
static constexpr double electron_Compton_length = hbarc / electron_mass_c2;
static constexpr double Bohr_radius =
    electron_Compton_length / fine_structure_const;

static constexpr double alpha_rcl2 =
    fine_structure_const * classic_electr_radius * classic_electr_radius;

static constexpr double twopi_mc2_rcl2 =
    twopi * electron_mass_c2 * classic_electr_radius * classic_electr_radius;
//
//
//
static constexpr double k_Boltzmann = 8.617385e-11 * MeV / kelvin;

//
//
//
static constexpr double STP_Temperature = 273.15 * kelvin;
static constexpr double STP_Pressure = 1. * atmosphere;
static constexpr double kGasThreshold = 10. * mg / cm3;

//
//
//
static constexpr double universe_mean_density = 1.e-25 * g / cm3;
}

}

#endif /* HEED_PHYSICAL_CONSTANTS_H */