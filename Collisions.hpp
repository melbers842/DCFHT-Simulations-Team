#pragma once
//
// Collisions.hpp
//
// Collision frequencies appearing in Eq. 13.22 - 13.26:
//
//   nu_i    ionization              (source term in continuity, Eq. 13.22)
//   nu_e    electron momentum loss  (drag in Eq. 13.23)
//   nu_in   ion-neutral scattering  (Eq. 13.24, 13.25)
//   nu_cex  charge exchange         (Eq. 13.24, 13.25)
//
// ===========================================================================
// IMPORTANT -- READ BEFORE PUBLISHING ANY RESULTS
//
// The cross sections below are order-of-magnitude placeholders good enough to
// get the solver running and debugged. They are NOT research-quality. Before
// any result goes in a paper, replace them with tabulated xenon cross sections
// from LXCat (lxcat.net) or the Biagi database, interpolated properly.
//
// Every constant marked PLACEHOLDER is a number someone has to justify.
// ===========================================================================
//
// DCFHT Simulation Team
//

#include "Constants.hpp"

#include <cmath>

namespace collisions {

// Mean thermal speed of a Maxwellian electron population, m/s.
// Te is in eV.
inline double electronThermalSpeed(double Te_eV) {
    using namespace constants;
    return std::sqrt(8.0 * kTe_J(Te_eV) / (kPi * kElectronMass));
}

// Electron cyclotron frequency, rad/s. B is the field magnitude in tesla.
inline double cyclotronFrequency(double B) {
    using namespace constants;
    return kElementaryCharge * B / kElectronMass;
}

// --- ionization -----------------------------------------------------------
//
// Rate coefficient <sigma v> for single ionization of xenon by electron
// impact, m^3/s, as a function of electron temperature in eV.
//
// PLACEHOLDER: simple Arrhenius-style fit. Correct order of magnitude across
// roughly 2-40 eV, but not accurate enough for publication.
inline double ionizationRateCoefficient(double Te_eV) {
    if (Te_eV <= 0.0) return 0.0;
    constexpr double kPrefactor = 1.8e-13;  // PLACEHOLDER, m^3/s
    return kPrefactor * std::sqrt(Te_eV)
         * std::exp(-constants::kXenonIonizationEnergy / Te_eV);
}

// Ionization frequency nu_i, 1/s. This is the nu_i multiplying n_e in the
// continuity source term of Eq. 13.22.
inline double ionizationFrequency(double n_n, double Te_eV) {
    return n_n * ionizationRateCoefficient(Te_eV);
}

// --- electron momentum transfer -------------------------------------------
//
// Classical electron-neutral momentum-transfer frequency, 1/s.
inline double electronNeutralFrequency(double n_n, double Te_eV) {
    constexpr double kSigmaEn = 2.7e-19;  // PLACEHOLDER, m^2
    return n_n * kSigmaEn * electronThermalSpeed(Te_eV);
}

// Anomalous ("Bohm") electron collision frequency, 1/s.
//
// Classical collisions underpredict cross-field electron transport in Hall
// thrusters by one to two orders of magnitude. Every working Hall thruster
// code adds an anomalous term; the usual first cut is a fixed fraction of the
// cyclotron frequency, alpha ~ 1/16. This is an open research problem in the
// field, not a shortcoming of your code -- but the value of alpha WILL change
// your answers, so treat it as a parameter to sweep, not a constant to forget.
inline double anomalousFrequency(double B, double alpha = 1.0 / 16.0) {
    return alpha * cyclotronFrequency(B);
}

// Total nu_e for the drag term in Eq. 13.23.
inline double electronMomentumFrequency(double n_n, double Te_eV, double B,
                                        double alpha = 1.0 / 16.0) {
    return electronNeutralFrequency(n_n, Te_eV) + anomalousFrequency(B, alpha);
}

// --- heavy-species collisions ---------------------------------------------
//
// Ion-neutral momentum transfer, 1/s. relSpeed is |v_i - v_n|.
inline double ionNeutralFrequency(double n_n, double relSpeed) {
    constexpr double kSigmaIn = 1.0e-18;  // PLACEHOLDER, m^2
    return n_n * kSigmaIn * relSpeed;
}

// Resonant charge exchange, Xe+ on Xe, 1/s.
inline double chargeExchangeFrequency(double n_n, double relSpeed) {
    constexpr double kSigmaCex = 5.0e-19;  // PLACEHOLDER, m^2
    return n_n * kSigmaCex * relSpeed;
}

}  // namespace collisions
