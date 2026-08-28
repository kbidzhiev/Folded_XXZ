#include "itensor/all.h"
#include <iostream>
#include <folded_xxz/initial_state.h>
#include <folded_xxz/output.h>
#include <folded_xxz/parameters.h>
#include <folded_xxz/profile.h>
#include <folded_xxz/simulation_config.h>
#include <folded_xxz/three_site_hamiltonian.h>
#include <folded_xxz/trotter_evolution.h>

void evolve_step(itensor::MPS &psi, int step, const SimulationSchedule &schedule, int center,
                 const itensor::Args &real_time_args, const itensor::Args &thermal_args,
                 TrotterEvolution &thermal_evolution, TrotterEvolution &half_thermal_evolution,
                 TrotterEvolution &real_time_evolution) {
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

int run_simulation(int argc, char *argv[]) {
  LOG_DURATION("MAIN");
  ThreeSiteParam param;
  param.parse_arguments(argc, argv); // Now param contains the parameters, default values or those
                                     // provided on the command-line

  param.print(std::cout); // Print parameters
  std::cout.precision(15);
  const auto config = make_simulation_config(param);
  const int N = config.site_count;

  itensor::SpinHalf sites(N, {"ConserveQNs=", false}); // HILBERT_SPACE = itensor::SpinHalf
  std::cout << "Constructing Hamiltonian" << std::endl;
  ThreeSiteHamiltonian hamiltonian(sites, param);
  const int dot = hamiltonian.center();
  auto H = hamiltonian.mpo();
  std::cout << "Finished constructing Hamiltonian" << std::endl;
  auto psi = create_initial_state(sites);
  std::cout << "N= " << N << std::endl;
  const auto energy = itensor::inner(psi, H, psi);
  std::cout << "2. Initial energy=" << energy << " .Norm = " << itensor::inner(psi, psi)
            << std::endl;

  auto args = itensor::Args("Method=", "DensityMatrix", "Cutoff", config.real_time_cutoff, "MaxDim",
                            config.max_bond_dimension, "Normalize",
                            false); // arguments for time dynamics

  auto args0 = itensor::Args("Method=", "DensityMatrix", "Cutoff", config.thermal_cutoff, "MaxDim",
                             config.max_bond_dimension, "Normalize",
                             false); // arguments for IMAGINARY time == initial temperature state

  auto output_files = open_output_files(config.output, dot);
  std::cout << "Trotter Gates for beta " << std::endl;
  const auto &schedule = config.schedule;
  TrotterEvolution expH_beta(sites, param, 0.5 * schedule.dbeta, EvolutionMode::ImaginaryTime);

  std::cout << "Trotter Gates Half for beta " << std::endl;
  TrotterEvolution expH_beta_half(sites, param, 0.5 * schedule.dbeta, EvolutionMode::ImaginaryTime,
                                  ThermalRegion::RightHalf);

  std::cout << "Trotter Gates for tau" << std::endl;
  param.set("hL", 0);
  param.set("hR", 0);
  TrotterEvolution expH(sites, param, itensor::Cplx_i * schedule.tau, EvolutionMode::RealTime);

  param.set("hL", 0);
  param.set("hR", 0);
  std::cout << "Finished preparing half-chain gates" << std::endl;

  for (int n = 0; n <= schedule.total_steps; ++n) {
    double time = (n - schedule.beta_steps_max) * schedule.tau;

    if (n < schedule.beta_steps_max) {
      time = n * schedule.dbeta;
      std::cout << "Beta(1/T) step #" << n << "/" << schedule.total_steps << "\t beta=" << time
                << std::endl;
    } else {
      std::cout << "Time step #" << n << "/" << schedule.total_steps << "\ttime=" << time
                << std::endl;
    }
    std::cout.flush();
    write_observables(psi, sites, H, config.output, output_files, n, time, schedule, N, dot);
    if (n < schedule.total_steps)
      evolve_step(psi, n, schedule, dot, args, args0, expH_beta, expH_beta_half, expH);
    std::cout << "max bond dim = " << itensor::maxLinkDim(psi) << std::endl;
    std::cout << "Norm = " << std::real(itensor::innerC(psi, psi)) << std::endl;
    std::cout << "Energy = " << std::real(itensor::innerC(psi, H, psi)) << std::endl << std::endl;
  }

  std::cout << "Final observables: \n"
            << "max bond dim = " << itensor::maxLinkDim(psi) << std::endl
            << "Energy = " << std::real(itensor::innerC(psi, H, psi)) << std::endl
            << std::endl;

  std::cout << "\nTime evolution complete.\n Done ! \n";

  return 0;
}
