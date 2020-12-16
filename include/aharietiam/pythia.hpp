#ifndef AHARIETIAM_PYTHIA_HPP
#define AHARIETIAM_PYTHIA_HPP

#include "aharietiam/constants.hpp"
#include "aharietiam/rambo.hpp"
#include "aharietiam/widths.hpp"
#include <Pythia8/Pythia.h>
#include <boost/histogram.hpp>
#include <boost/histogram/axis/regular.hpp>
#include <boost/histogram/weight.hpp>
#include <string>

using Pythia8::Pythia;
using Pythia8::Vec4;
using namespace boost::histogram;

namespace aharietiam {

/**
 * Class for generating gamma-ray and neutrino spectra from the decay
 * of a RH-Neutrino.
 */
class PythiaDecaySpectum {
private:
  typedef axis::regular<double, axis::transform::log> log_axis;
  typedef histogram<std::tuple<log_axis>> hist;

  const double m_mass;
  const double m_stheta;

  const size_t m_nevents;
  const size_t m_nbins;
  const double m_emin;
  const double m_emax;

  double bf_h_nu;
  double bf_z_nu;
  double bf_w_ell;

  int lep_id;
  double lep_mass;

  hist photons;
  hist nu_e;
  hist nu_mu;
  hist nu_tau;

  Pythia pythia{};

public:
  PythiaDecaySpectum(double mass, double stheta, const std::string &lepton,
                     size_t nevents, size_t nbins, double emin, double emax)
      : m_mass(mass), m_stheta(stheta), m_nevents(nevents), m_nbins(nbins),
        m_emin(emin), m_emax(emax) {

    // Key requirement: switch off ProcessLevel, and thereby also PartonLevel.
    pythia.readString("ProcessLevel:all = off");

    // Switch off automatic event listing in favour of manual.
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");

    pythia.readString("13:mayDecay = on");  // Turn on charged muon decay
    pythia.readString("130:mayDecay = on"); // Turn on long kaon decay
    pythia.readString("211:mayDecay = on"); // Turn on charged pion decay
    pythia.readString("321:mayDecay = on"); // Turn on charged kaon decay
    pythia.init();

    if (lepton == "e") {
      lep_id = kELECTRON_ID;
      lep_mass = kELECTRON_MASS;
    } else if (lepton == "mu") {
      lep_id = kMUON_ID;
      lep_mass = kMUON_MASS;
    } else {
      lep_id = kTAU_ID;
      lep_mass = kTAU_MASS;
    }

    const double bf_h_nu = branching_fraction_h_nu(mass, stheta, lep_mass);
    const double bf_z_nu = branching_fraction_z_nu(mass, stheta, lep_mass);
    const double bf_w_ell = branching_fraction_w_ell(mass, stheta, lep_mass);

    photons = make_histogram(axis::regular<double, axis::transform::log>{
        static_cast<unsigned int>(nbins), emin, emax});
    nu_e = make_histogram(axis::regular<double, axis::transform::log>{
        static_cast<unsigned int>(nbins), emin, emax});
    nu_mu = make_histogram(axis::regular<double, axis::transform::log>{
        static_cast<unsigned int>(nbins), emin, emax});
    nu_tau = make_histogram(axis::regular<double, axis::transform::log>{
        static_cast<unsigned int>(nbins), emin, emax});
  }

private:
  auto fill_particles(std::vector<int> &ids, std::vector<Vec4> &momenta)
      -> void;
  auto simulate_h_nu() -> void;
  auto simulate_z_nu() -> void;
  auto simulate_wp_lm() -> void;
  auto simulate_wm_lp() -> void;
};

/**
 * Add particles into the Pythia event record.
 * @param ids
 * @param momenta
 * @param event
 * @param rndm
 */
void PythiaDecaySpectum::fill_particles(std::vector<int> &ids,
                                        std::vector<Vec4> &momenta) {
  for (size_t i = 0; i < ids.size(); i++) {
    int id = ids[i];
    int iNew;
    // treat colored states
    if (id == 1 || id == 2 || id == 3 || id == 4 || id == 5 || id == 6) {
      iNew = pythia.event.append(id, 23, 101, 0, momenta[i].px(),
                                 momenta[i].py(), momenta[i].pz(),
                                 momenta[i].e(), momenta[i].mCalc());
    } else if (id == -1 || id == -2 || id == -3 || id == -4 || id == -5 ||
               id == -6) {
      iNew = pythia.event.append(id, 23, 0, 101, momenta[i].px(),
                                 momenta[i].py(), momenta[i].pz(),
                                 momenta[i].e(), momenta[i].mCalc());
    } else {
      iNew = pythia.event.append(id, 23, 0, 0, momenta[i].px(), momenta[i].py(),
                                 momenta[i].pz(), momenta[i].e(),
                                 momenta[i].mCalc());
    }

    // Generate a lifetime, to give decay away from a primary vertex.
    if (id != 11 && id != -11 && id != 12 && id != 14 && id != 15 && id != 22)
      pythia.event[iNew].tau(pythia.event[iNew].tau0() * pythia.rndm.exp());
  }
}

auto PythiaDecaySpectum::simulate_h_nu() -> void { pythia.event.reset(); }

void run_event(Pythia &pythia, PythiaRun &pythia_run, const double mass,
               const double m1, const double m2, const double bf, const int id1,
               const int id2) {
  // Reset event record to allow for a new event.
  pythia.event.reset();

  std::std::vector<int> fsp_ids = {id1, id2};
  double p = sqrt((mass - m1 - m2) * (mass + m1 - m2) * (mass - m1 + m2) *
                  (mass + m1 + m2)) /
             (2.0 * mass);
  Vec4 p1{sqrt(p * p + m1 * m1), 0.0, 0.0, p};
  Vec4 p2{sqrt(p * p + m2 * m2), 0.0, 0.0, -p};
  std::std::vector<Vec4> momenta = {p1, p2};
  fill_particles(fsp_ids, momenta, pythia.event, pythia.rndm);

  // Generate events. Quit if failure.
  if (!pythia.next()) {
    cout << " Event generation aborted prematurely, owing to error!\n";
    return;
  }

  // Loop over all particles.
  for (int i = 0; i < pythia.event.size(); ++i) {
    if (pythia.event[i].isFinal()) {
      int idAbs = pythia.event[i].idAbs();
      double eI = pythia.event[i].e();
      double _weight = pythia.info.weight() * bf;
      if (idAbs == 22) {
        pythia_run.photons(eI, weight(_weight));
      } else if (idAbs == 12) {
        pythia_run.nu_e(eI, weight(_weight));
      } else if (idAbs == 14) {
        pythia_run.nu_mu(eI, weight(_weight));
      } else if (idAbs == 16) {
        pythia_run.nu_tau(eI, weight(_weight));
      }
    }
  }
}

/**
 * Generate a the spectra into photon and neutrinos from the decay of a
 * RH-Neutrino.
 * @param pythia_run
 * @param mass
 * @param stheta
 * @param ml
 */
void generate_spectrum(PythiaRun &pythia_run, const double mass,
                       const double stheta, const std::string &lepton) {

  // Initialize pythia
  Pythia pythia{};
  // Key requirement: switch off ProcessLevel, and thereby also PartonLevel.
  pythia.readString("ProcessLevel:all = off");
  // Switch off automatic event listing in favour of manual.
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  // Turn on charged muon decay
  pythia.readString("13:mayDecay = on");
  // Turn on long kaon decay
  pythia.readString("130:mayDecay = on");
  // Turn on charged pion decay
  pythia.readString("211:mayDecay = on");
  // Turn on charged kaon decay
  pythia.readString("321:mayDecay = on");
  pythia.init();

  int lep_id = kELECTRON_ID;
  int nu_id = kELECTRON_NEUTRINO_ID;
  double ml = kELECTRON_MASS;
  if (lepton == "mu") {
    lep_id = kMUON_ID;
    nu_id = kMUON_NEUTRINO_ID;
    ml = kMUON_MASS;
  } else {
    lep_id = kTAU_ID;
    nu_id = kTAU_NEUTRINO_ID;
    ml = kTAU_MASS;
  }

  const double bf_h_nu = branching_fraction_h_nu(mass, stheta, ml);
  const double bf_z_nu = branching_fraction_z_nu(mass, stheta, ml);
  const double bf_w_ell = branching_fraction_w_ell(mass, stheta, ml);

  for (size_t iEvent = 0; iEvent < pythia_run.num_events; iEvent += 1) {
    run_event(pythia, pythia_run, mass, kHIGGS_MASS, 0.0, bf_h_nu, kHIGGS_ID,
              nu_id);
    run_event(pythia, pythia_run, mass, kZ_BOSON_MASS, 0.0, bf_z_nu,
              kZ_BOSON_ID, nu_id);
    run_event(pythia, pythia_run, mass, kW_BOSON_MASS, ml, bf_w_ell,
              -kW_BOSON_ID, lep_id);
    run_event(pythia, pythia_run, mass, kW_BOSON_MASS, ml, bf_w_ell,
              kW_BOSON_ID, -lep_id);
  }

  // Normalize the events
  for (int i = 0; i < pythia_run.photons.axis().size(); i++) {
    pythia_run.photons.at(i) /= double(pythia_run.num_events);
  }
  for (int i = 0; i < pythia_run.nu_e.axis().size(); i++) {
    pythia_run.nu_e.at(i) /= double(pythia_run.num_events);
  }
  for (int i = 0; i < pythia_run.nu_mu.axis().size(); i++) {
    pythia_run.nu_mu.at(i) /= double(pythia_run.num_events);
  }
  for (int i = 0; i < pythia_run.nu_tau.axis().size(); i++) {
    pythia_run.nu_tau.at(i) /= double(pythia_run.num_events);
  }
}

} // namespace aharietiam

#endif // AHARIETIAM_PYTHIA_HPP
