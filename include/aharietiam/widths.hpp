/**
 * widths.hpp
 * This file contains definitions for the decay widths of the RH-Neutrino into
 * 2-body final states.
 */

#ifndef AHARIETIAM_WIDTHS_HPP
#define AHARIETIAM_WIDTHS_HPP

#include "aharietiam/constants.hpp"
#include <cmath>

namespace aharietiam {

/**
 * Compute the partial width of a RH neutrino decaying into a Higgs and
 * neutrino.
 * @param mass Mass of the RH neutrino.
 * @param stheta Mixing angle between left- and right-handed neutrinos.
 * @param ml Mass of the lepton associate with the active neutrino.
 */
auto width_h_nu(double mass, double stheta, double /*ml*/) -> double {
  if (mass > kHIGGS_MASS) {
    const double MH2 = kHIGGS_MASS * kHIGGS_MASS;
    const double MN2 = mass * mass;
    return ((MH2 - MN2) * sqrt(pow(MH2 - MN2, 2) / MN2) *
            pow(-stheta + 2 * pow(stheta, 3), 2)) /
           (16.0 * M_PI * (-1 + stheta * stheta) * kHIGGS_VEV * kHIGGS_VEV);
  }
  return 0.0;
}

/**
 * Compute the partial width of a RH neutrino decaying into a Z-Boson and
 * neutrino.
 * @param mass Mass of the RH neutrino.
 * @param stheta Mixing angle between left- and right-handed neutrinos.
 * @param ml Mass of the lepton associate with the active neutrino.
 */
auto width_z_nu(double mass, double stheta, double /*ml*/) -> double {
  if (mass > kZ_BOSON_MASS) {
    constexpr double MZ2 = kZ_BOSON_MASS * kZ_BOSON_MASS;
    constexpr double MZ4 = MZ2 * MZ2;
    constexpr double CW2 = kCOS_THETA_WEAK * kCOS_THETA_WEAK;
    constexpr double SW2 = kSIN_THETA_WEAK_SQRD;
    const double MN2 = mass * mass;
    const double MN4 = MN2 * MN2;
    return -(kALPHA_EM * fabs(MN2 - MZ2) / mass * (MN4 + MN2 * MZ2 - 2 * MZ4) *
             (-1.0 + stheta) * stheta * stheta * (1.0 + stheta)) /
           (16.0 * CW2 * MN2 * MZ2 * SW2);
  }
  return 0.0;
}

/**
 * Compute the partial width of a RH neutrino decaying into a W-Boson and
 * lepton.
 * @param mass Mass of the RH neutrino.
 * @param stheta Mixing angle between left- and right-handed neutrinos.
 * @param ml Mass of the lepton associate with the active neutrino.
 */
auto width_w_ell(double mass, double stheta, double ml) -> double {
  if (mass > kW_BOSON_MASS + ml) {
    constexpr double MZ2 = kZ_BOSON_MASS * kZ_BOSON_MASS;
    constexpr double MZ4 = MZ2 * MZ2;
    constexpr double CW2 = kCOS_THETA_WEAK * kCOS_THETA_WEAK;
    constexpr double SW2 = kSIN_THETA_WEAK_SQRD;
    const double MN2 = mass * mass;
    const double MN4 = MN2 * MN2;
    return -(kALPHA_EM * fabs(MN2 - MZ2) / mass * (MN4 + MN2 * MZ2 - 2 * MZ4) *
             (-1.0 + stheta) * stheta * stheta * (1.0 + stheta)) /
           (16.0 * CW2 * MN2 * MZ2 * SW2);
  }
  return 0.0;
}

/**
 * Compute the total width of a RH neutrino.
 * @param mass Mass of the RH neutrino.
 * @param stheta Mixing angle between left- and right-handed neutrinos.
 * @param ml Mass of the lepton associate with the active neutrino.
 * @return width The total width.
 */
auto width_tot(double mass, double stheta, double ml) -> double {
  return width_h_nu(mass, stheta, ml) + width_z_nu(mass, stheta, ml) +
         width_w_ell(mass, stheta, ml);
}

/**
 * Compute the branching fraction of a RH neutrino decaying into a Higgs and
 * neutrino.
 * @param mass Mass of the RH neutrino.
 * @param stheta Mixing angle between left- and right-handed neutrinos.
 * @param ml Mass of the lepton associate with the active neutrino.
 */
auto branching_fraction_h_nu(double mass, double stheta, double ml) -> double {
  double w = width_h_nu(mass, stheta, ml);
  if (w > 0.0) {
    double width =
        w + width_z_nu(mass, stheta, ml) + width_w_ell(mass, stheta, ml);
    return w / width;
  }
  return 0.0;
}

/**
 * Compute the branching fraction of a RH neutrino decaying into a Z-Boson and
 * neutrino.
 * @param mass Mass of the RH neutrino.
 * @param stheta Mixing angle between left- and right-handed neutrinos.
 * @param ml Mass of the lepton associate with the active neutrino.
 */
auto branching_fraction_z_nu(double mass, double stheta, double ml) -> double {
  double w = width_z_nu(mass, stheta, ml);
  if (w > 0.0) {
    double width =
        w + width_h_nu(mass, stheta, ml) + width_w_ell(mass, stheta, ml);
    return w / width;
  }
  return 0.0;
}

/**
 * Compute the branching fraction of a RH neutrino decaying into a W-Boson and
 * lepton.
 * @param mass Mass of the RH neutrino.
 * @param stheta Mixing angle between left- and right-handed neutrinos.
 * @param ml Mass of the lepton associate with the active neutrino.
 */
auto branching_fraction_w_ell(double mass, double stheta, double ml) -> double {
  double w = width_w_ell(mass, stheta, ml);
  if (w > 0.0) {
    double width =
        w + width_h_nu(mass, stheta, ml) + width_z_nu(mass, stheta, ml);
    return w / width;
  }
  return 0.0;
}

} // namespace aharietiam
#endif // AHARIETIAM_WIDTHS_HPP
