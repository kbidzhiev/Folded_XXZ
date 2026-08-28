#include "itensor/all.h"

#include <folded_xxz/model_config.h>
#include <folded_xxz/trotter_evolution.h>

#include <cmath>
#include <iostream>

namespace {

void expect(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << std::endl;
    std::exit(1);
  }
}

} // namespace

int main() {
  itensor::SpinHalf sites(12, {"ConserveQNs=", false});
  const ModelConfig config{12, 1.0, 0.0, 0.0, 2};

  const TrotterEvolution full_thermal(sites, config, 0.01, EvolutionMode::ImaginaryTime);
  const TrotterEvolution right_half_thermal(sites, config, 0.01, EvolutionMode::ImaginaryTime,
                                            ThermalRegion::RightHalf);
  expect(full_thermal.gate_count() > right_half_thermal.gate_count(),
         "right-half thermal evolution has fewer gates");

  itensor::InitState initial_state(sites, "Up");
  itensor::MPS psi(initial_state);
  const TrotterEvolution real_time(sites, config, itensor::Cplx_i * 0.01, EvolutionMode::RealTime);
  real_time.evolve(psi, {"Cutoff", 1e-12, "MaxDim", 100, "Normalize", false});
  expect(std::abs(itensor::innerC(psi, psi).real() - 1.0) < 1e-10,
         "real-time evolution preserves normalization");
}
