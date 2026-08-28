#include "itensor/all.h"

#include <folded_xxz/parameters.h>
#include <folded_xxz/trotter_evolution.h>

#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>

namespace {

void expect(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << std::endl;
    std::exit(1);
  }
}

void expect_invalid_argument(const std::function<void()> &operation, const char *message) {
  try {
    operation();
  } catch (const std::invalid_argument &) {
    return;
  }
  expect(false, message);
}

} // namespace

int main() {
  itensor::SpinHalf sites(12, {"ConserveQNs=", false});
  ThreeSiteParam parameters;

  const TrotterEvolution full_thermal(sites, parameters, 0.01, EvolutionMode::ImaginaryTime);
  const TrotterEvolution right_half_thermal(sites, parameters, 0.01, EvolutionMode::ImaginaryTime,
                                            ThermalRegion::RightHalf);
  expect(full_thermal.gate_count() > right_half_thermal.gate_count(),
         "right-half thermal evolution has fewer gates");

  itensor::InitState initial_state(sites, "Up");
  itensor::MPS psi(initial_state);
  const TrotterEvolution real_time(sites, parameters, itensor::Cplx_i * 0.01,
                                   EvolutionMode::RealTime);
  real_time.evolve(psi, {"Cutoff", 1e-12, "MaxDim", 100, "Normalize", false});
  expect(std::abs(itensor::innerC(psi, psi).real() - 1.0) < 1e-10,
         "real-time evolution preserves normalization");

  parameters.set("TrotterOrder", 3);
  expect_invalid_argument(
      [&] { TrotterEvolution(sites, parameters, 0.01, EvolutionMode::ImaginaryTime); },
      "invalid Trotter order");
}
