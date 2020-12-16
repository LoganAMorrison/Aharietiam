#include <boost/format.hpp>
#include <gtest/gtest.h>
#include <rhndecay/pythia.hpp>
#include <string>

using namespace rhndecay;

TEST(TestPythia, TestHighEnergy) {
  double mass = 1e3;
  double stheta = 1e-3;
  string lepton = "e";

  PythiaRun run{5000, 15, 1.0, mass / 2.0};

  generate_spectrum(run, mass, stheta, lepton);
  std::ostringstream os;
  for (auto &&x : indexed(run.photons)) {
    const auto i = x.index(0); // current index along first axis
    const auto j = x.index(1); // current index along second axis
    const auto b0 = x.bin(0);  // current bin interval along first axis
    const auto b1 = x.bin(1);  // current bin interval along second axis
    const auto v = *x;         // "dereference" to get the bin value
    os << boost::format("%i %i [%2i, %i) [%2i, %i): %i\n") % i % j %
              b0.lower() % b0.upper() % b1.lower() % b1.upper() % v;
  }
}
