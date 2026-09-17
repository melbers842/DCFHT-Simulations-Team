#pragma once
//
// Constants.hpp
//
// Physical constants and xenon propellant properties, in SI units throughout.
//
// Unit convention for the whole project: SI everywhere, EXCEPT electron
// temperature, which is stored in electron-volts because every formula in the
// plasma literature is written that way. Convert with kTe_J() when a formula
// needs joules. Getting this wrong is the single most common source of
// silently wrong answers in plasma codes -- keep the eV/J boundary at this
// one place.
//
// DCFHT Simulation Team
//

namespace constants {

// --- universal ------------------------------------------------------------
inline constexpr double kElementaryCharge = 1.602176634e-19;  // C
inline constexpr double kBoltzmann        = 1.380649e-23;     // J/K
inline constexpr double kElectronMass     = 9.1093837015e-31; // kg
inline constexpr double kEpsilon0         = 8.8541878128e-12; // F/m
inline constexpr double kAmu              = 1.66053906660e-27;// kg
inline constexpr double kPi               = 3.14159265358979323846;

// --- xenon propellant -----------------------------------------------------
inline constexpr double kXenonMass        = 131.293 * kAmu;   // kg
inline constexpr double kXenonIonizationEnergy = 12.1298;     // eV (E_i in Eq. 13.26)

// --- lecture-specified conditions -----------------------------------------
// Neutral temperature is held constant, per Section 13.4.
inline constexpr double kNeutralTemperature = 300.0;          // K

// --- conversions ----------------------------------------------------------
// Electron temperature is carried in eV. This turns it into joules, i.e. k*T_e.
inline constexpr double kTe_J(double Te_eV) {
    return Te_eV * kElementaryCharge;
}

// k*T_n in joules, for the neutral pressure term in Eq. 13.25.
inline constexpr double kTn_J() {
    return kBoltzmann * kNeutralTemperature;
}

}  // namespace constants
