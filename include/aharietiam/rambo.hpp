//
// Created by Logan Morrison on 12/13/20.
//

#ifndef AHARIETIAM_RAMBO_HPP
#define AHARIETIAM_RAMBO_HPP

#include "Pythia8/Pythia.h"
#include <cmath>
#include <functional>
#include <mutex>
#include <numeric>
#include <random>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

using Pythia8::Vec4;

namespace aharietiam {

struct PhaseSpaceEvent {
  std::vector<Vec4> momenta;
  double weight;
};

/**
 * Uniform random number generator that is thread-safe.
 * @return random number between (0,1)
 */
static auto phase_space_uniform_rand() -> double {
  static thread_local std::random_device rd{};
  static thread_local std::mt19937 generator{rd()};
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  return distribution(generator);
}

/* Function for finding the scaling parameter to turn mass-less four-vectors
 * into four-vectors with the correct masses.
 * @param momenta 4-momenta of final-state particles
 */
auto compute_scale_factor(const std::vector<Vec4> &momenta,
                          const std::vector<double> &fsp_masses,
                          const double cme) -> double {
  static thread_local const int MAX_ITER = 50;
  static thread_local const double TOL = 1e-4;

  double mass_sum = std::accumulate(fsp_masses.begin(), fsp_masses.end(), 0.0);

  double xi = sqrt(1.0 - (mass_sum / cme) * (mass_sum / cme));

  int iter_count = 0;
  bool converged = false;
  do { // Perform newton iterations to solve for xi
    double f = -cme;
    double df = 0.0;

    for (size_t i = 0; i < fsp_masses.size(); i++) {
      // Compute residual and derivative of residual
      double m2 = fsp_masses[i] * fsp_masses[i];
      double xi2 = xi * xi;
      double e2 = momenta[i].e() * momenta[i].e();
      double del_f = sqrt(m2 + xi2 * e2);
      f += del_f;
      df += xi * e2 / del_f;
    }

    // Newton correction
    double delta_xi = -(f / df);
    xi += delta_xi;

    iter_count++;
    if (fabs(delta_xi) < TOL || iter_count >= MAX_ITER) {
      converged = true;
    }
  } while (!converged);
  return xi;
}

/**
 * Initialize the four-momenta with isotropic, random four-momenta with
 * energies, q₀, distributed according to q₀ * exp(-q₀).
 * @param momenta 4-momenta of final-state particles
 */
auto initialize_four_momenta(std::vector<Vec4> *momenta) -> void {

  std::uniform_real_distribution<double> distribution(0.0, 1.0);

  for (auto &mom : *momenta) {
    double rho1 = phase_space_uniform_rand();
    double rho2 = phase_space_uniform_rand();
    double rho3 = phase_space_uniform_rand();
    double rho4 = phase_space_uniform_rand();

    double c = 2.0 * rho1 - 1.0;
    double phi = 2.0 * M_PI * rho2;

    mom.e(-log(rho3 * rho4));
    mom.px(mom.e() * sqrt(1.0 - c * c) * cos(phi));
    mom.py(mom.e() * sqrt(1.0 - c * c) * sin(phi));
    mom.pz(mom.e() * c);
  }
}

/**
 * Boost the four-momenta into the center-of-mass frame and compute the
 * initial weight of the event.
 * @param momenta 4-momenta of final-state particles
 */
auto boost_four_momenta(std::vector<Vec4> *momenta, const double cme) -> void {
  // Total momentum and its mass
  Vec4 Q = std::accumulate(momenta->begin(), momenta->end(), Vec4{});
  double massQ = Q.mCalc();

  // Boost three-vector
  double bx = -Q.px() / massQ;
  double by = -Q.py() / massQ;
  double bz = -Q.pz() / massQ;
  // Boost factors
  double x = cme / massQ;
  double gamma = Q.e() / massQ;
  double a = 1.0 / (1.0 + gamma);

  for (auto &mom : *momenta) {
    double qe = mom.e();
    double qx = mom.px();
    double qy = mom.py();
    double qz = mom.pz();

    double b_dot_q = bx * qx + by * qy + bz * qz;

    mom.e(x * (gamma * qe + b_dot_q));
    mom.px(x * (qx + bx * qe + a * b_dot_q * bx));
    mom.py(x * (qy + by * qe + a * b_dot_q * by));
    mom.pz(x * (qz + bz * qe + a * b_dot_q * bz));
  }
}

/**
 * Correct the masses of the four-momenta and correct the weight of the
 * event.
 * @param momenta 4-momenta of final-state particles
 * @return new event weight factor
 */
auto correct_masses(std::vector<Vec4> *momenta,
                    const std::vector<double> &fsp_masses, const double cme)
    -> double {
  double xi = compute_scale_factor(*momenta, fsp_masses, cme);

  double term1 = 0.0;
  double term2 = 0.0;
  double term3 = 1.0;

  for (size_t i = 0; i < fsp_masses.size(); i++) {
    double m = fsp_masses[i];
    double eng = momenta->at(i).e();
    momenta->at(i).e(sqrt(m * m + (xi * eng) * (xi * eng)));
    momenta->at(i).px(momenta->at(i).px() * xi);
    momenta->at(i).py(momenta->at(i).py() * xi);
    momenta->at(i).pz(momenta->at(i).pz() * xi);

    double mod = sqrt(momenta->at(i).px() * momenta->at(i).px() +
                      momenta->at(i).py() * momenta->at(i).py() +
                      momenta->at(i).pz() * momenta->at(i).pz());
    eng = momenta->at(i).e();

    term1 += mod / cme;
    term2 += mod * mod / eng;
    term3 *= mod / eng;
  }

  term1 = pow(term1, 2.0 * fsp_masses.size() - 3.0);
  term2 = 1.0 / term2;

  // re-weight
  return term1 * term2 * term3 * cme;
}

/**
 * generate single phase space event
 * @return event
 */
auto generate_event(const std::vector<double> &fsp_masses, const double cme)
    -> PhaseSpaceEvent {
  std::vector<Vec4> momenta(fsp_masses.size(), Vec4{});
  const auto num_fsp_d = static_cast<double>(fsp_masses.size());
  const double base_weight = pow(M_PI / 2.0, num_fsp_d - 1.0) *
                             pow(cme, 2.0 * num_fsp_d - 4.0) /
                             tgamma(num_fsp_d) / tgamma(num_fsp_d - 1.0) *
                             pow(2.0 * M_PI, 4.0 - 3.0 * num_fsp_d);
  initialize_four_momenta(&momenta);
  boost_four_momenta(&momenta, cme);
  const double weight = correct_masses(&momenta, fsp_masses, cme) * base_weight;

  return PhaseSpaceEvent{std::move(momenta), weight};
}

} // namespace aharietiam

#endif // AHARIETIAM_RAMBO_HPP
