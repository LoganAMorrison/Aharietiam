#ifndef AHARIETIAM_CONSTANTS_HPP
#define AHARIETIAM_CONSTANTS_HPP

namespace aharietiam {

// ==========================
// ---- Physical Constants --
// ==========================

// Various physical constants
static constexpr double kG_FERMI = 1.1663787e-5;
static constexpr double kHIGGS_VEV = 246.21965;
static constexpr double kALPHA_EM = 1.0 / 137.0; // at p^2 = 0
static constexpr double kSIN_THETA_WEAK = 0.480853;
static constexpr double kSIN_THETA_WEAK_SQRD = 0.23122;
static constexpr double kCOS_THETA_WEAK = 0.876801;
static constexpr double kM_PLANK = 1.220910e19;
static constexpr double kRHO_CRIT = 1.05375e-5;
static constexpr double kS_TODAY = 2891.2;
static constexpr double kT_CMB = 2.56215e-10;
static constexpr double kT_BBN = 0.0001; // 0.1 MeV in GeV
static constexpr double kOMEGA_H2_CDM = 0.1198;

// Masses
static constexpr double kELECTRON_MASS = 0.5109989461e-3;
static constexpr double kMUON_MASS = 105.6583745e-3;
static constexpr double kTAU_MASS = 1776.86e-3;
static constexpr double kUP_QUARK_MASS = 2.16e-3;
static constexpr double kDOWN_QUARK_MASS = 4.67e-3;
static constexpr double kSTRANGE_QUARK_MASS = 93e-3;
static constexpr double kCHARM_QUARK_MASS = 1.27;
static constexpr double kBOTTOM_QUARK_MASS = 4.18;
static constexpr double kTOP_QUARK_MASS = 172.9;
static constexpr double kW_BOSON_MASS = 80.379;
static constexpr double kZ_BOSON_MASS = 91.1876;
static constexpr double kHIGGS_MASS = 125.10;
static constexpr double kNEUTRAL_PION_MASS = 134.9766e-3;
static constexpr double kCHARGED_PION_MASS = 139.57018e-3;
static constexpr double kNEUTRAL_KAON_MASS = 497.61e-3;
static constexpr double kLONG_KAON_MASS = 497.614e-3;
static constexpr double kCHARGED_KAON_MASS = 493.68e-3;

// Boson widths
static constexpr double kW_BOSON_WIDTH = 2.085;
static constexpr double kZ_BOSON_WIDTH = 2.4952;
static constexpr double kHIGGS_WIDTH = 4.07e-3;

// PDG Codes
static constexpr int kELECTRON_ID = 11;
static constexpr int kELECTRON_NEUTRINO_ID = 12;
static constexpr int kMUON_ID = 13;
static constexpr int kMUON_NEUTRINO_ID = 14;
static constexpr int kTAU_ID = 15;
static constexpr int kTAU_NEUTRINO_ID = 16;
static constexpr int kHIGGS_ID = 25;
static constexpr int kZ_BOSON_ID = 23;
static constexpr int kW_BOSON_ID = 24;
static constexpr int kUP_QUARK_ID = 2;
static constexpr int kCHARM_QUARK_ID = 4;
static constexpr int kTOP_QUARK_ID = 6;
static constexpr int kDOWN_QUARK_ID = 1;
static constexpr int kSTRANGE_QUARK_ID = 3;
static constexpr int kBOTTOM_QUARK_ID = 5;

double id_to_mass(const int id) {
  switch (id) {
  case kDOWN_QUARK_ID:
    return kDOWN_QUARK_MASS;
  case kUP_QUARK_ID:
    return kUP_QUARK_MASS;
  case kSTRANGE_QUARK_ID:
    return kSTRANGE_QUARK_MASS;
  case kCHARM_QUARK_ID:
    return kCHARM_QUARK_MASS;
  case kBOTTOM_QUARK_ID:
    return kBOTTOM_QUARK_MASS;
  case kTOP_QUARK_ID:
    return kTOP_QUARK_MASS;
  case kELECTRON_ID:
    return kELECTRON_MASS;
  case kMUON_ID:
    return kMUON_MASS;
  case kTAU_ID:
    return kTAU_MASS;
  case kW_BOSON_ID:
    return kW_BOSON_MASS;
  case kZ_BOSON_ID:
    return kZ_BOSON_MASS;
  case kHIGGS_ID:
    return kHIGGS_MASS;
  default:
    return 0.0;
  }
}

} // namespace aharietiam

#endif // AHARIETIAM_CONSTANTS_HPP
