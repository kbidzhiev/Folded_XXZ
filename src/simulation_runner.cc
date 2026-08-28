#include "itensor/all.h"
#include "output.h"
#include "profile.h"
#include "simulation_runner.h"

#include <iostream>
#include <utility>

#include <folded_xxz/initial_state.h>
#include <folded_xxz/three_site_hamiltonian.h>
#include <folded_xxz/trotter_evolution.h>

namespace {

void evolve_step(itensor::MPS &psi, int step, const SimulationSchedule &schedule, int center,
                 const itensor::Args &real_time_args, const itensor::Args &thermal_args,
                 const TrotterEvolution &thermal_evolution,
                 const TrotterEvolution &half_thermal_evolution,
                 const TrotterEvolution &real_time_evolution) {
  if (step < schedule.beta_steps_min) {
    std::cout << "Temperature evolution with H" << std::endl;
    thermal_evolution.evolve(psi, thermal_args);
    psi.orthogonalize(real_time_args);
    std::cout << "dot = " << center + 1 << std::endl;
    psi.normalize();
    return;
  }

  if (step < schedule.beta_steps_max) {
    std::cout << "Temperature evolution with half-chain H" << std::endl;
    std::cout << "n/beta_steps_max = " << step << "/" << schedule.beta_steps_max << std::endl;
    half_thermal_evolution.evolve(psi, thermal_args);
    psi.orthogonalize(real_time_args);
    psi.normalize();
    return;
  }

  std::cout << "Time evolution" << std::endl;
  real_time_evolution.evolve(psi, real_time_args);
  psi.orthogonalize(real_time_args);
}

} // namespace

SimulationRunner::SimulationRunner(SimulationConfig config) : config_(std::move(config)) {}

int SimulationRunner::run() {
  LOG_DURATION("MAIN");
  const int site_count = config_.model.site_count;
  itensor::SpinHalf sites(site_count, {"ConserveQNs=", false});
  std::cout << "Constructing Hamiltonian" << std::endl;
  ThreeSiteHamiltonian hamiltonian(sites, config_.model);
  const int center = hamiltonian.center();
  const auto mpo = hamiltonian.mpo();
  std::cout << "Finished constructing Hamiltonian" << std::endl;

  auto psi = create_initial_state(sites);
  std::cout << "N= " << site_count << std::endl;
  std::cout << "2. Initial energy=" << itensor::inner(psi, mpo, psi)
            << " .Norm = " << itensor::inner(psi, psi) << std::endl;

  const auto real_time_args =
      itensor::Args("Method=", "DensityMatrix", "Cutoff", config_.real_time_cutoff, "MaxDim",
                    config_.max_bond_dimension, "Normalize", false);
  const auto thermal_args =
      itensor::Args("Method=", "DensityMatrix", "Cutoff", config_.thermal_cutoff, "MaxDim",
                    config_.max_bond_dimension, "Normalize", false);

  const auto &schedule = config_.schedule;
  ObservableWriter output(config_.output, schedule, center);
  std::cout << "Trotter Gates for beta " << std::endl;
  const TrotterEvolution thermal_evolution(sites, config_.model, 0.5 * schedule.dbeta,
                                           EvolutionMode::ImaginaryTime);
  std::cout << "Trotter Gates Half for beta " << std::endl;
  const TrotterEvolution half_thermal_evolution(sites, config_.model, 0.5 * schedule.dbeta,
                                                EvolutionMode::ImaginaryTime,
                                                ThermalRegion::RightHalf);
  std::cout << "Trotter Gates for tau" << std::endl;
  const TrotterEvolution real_time_evolution(sites, config_.model, itensor::Cplx_i * schedule.tau,
                                             EvolutionMode::RealTime);
  std::cout << "Finished preparing half-chain gates" << std::endl;

  for (int step = 0; step <= schedule.total_steps; ++step) {
    double time = (step - schedule.beta_steps_max) * schedule.tau;
    if (step < schedule.beta_steps_max) {
      time = step * schedule.dbeta;
      std::cout << "Beta(1/T) step #" << step << "/" << schedule.total_steps << "\t beta=" << time
                << std::endl;
    } else {
      std::cout << "Time step #" << step << "/" << schedule.total_steps << "\ttime=" << time
                << std::endl;
    }
    std::cout.flush();
    output.write(psi, sites, mpo, step, time);
    if (step < schedule.total_steps) {
      evolve_step(psi, step, schedule, center, real_time_args, thermal_args, thermal_evolution,
                  half_thermal_evolution, real_time_evolution);
    }
    std::cout << "max bond dim = " << itensor::maxLinkDim(psi) << std::endl;
    std::cout << "Norm = " << std::real(itensor::innerC(psi, psi)) << std::endl;
    std::cout << "Energy = " << std::real(itensor::innerC(psi, mpo, psi)) << std::endl << std::endl;
  }

  std::cout << "Final observables: \n"
            << "max bond dim = " << itensor::maxLinkDim(psi) << std::endl
            << "Energy = " << std::real(itensor::innerC(psi, mpo, psi)) << std::endl
            << std::endl;
  std::cout << "\nTime evolution complete.\n Done ! \n";
  return 0;
}
