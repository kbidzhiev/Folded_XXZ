#include <folded_xxz/parameters.h>
#include <folded_xxz/simulation_config.h>
#include <folded_xxz/simulation_runner.h>

#include <iostream>

int run_simulation(int argc, char *argv[]) {
  ThreeSiteParam param;
  param.parse_arguments(argc, argv); // Now param contains the parameters, default values or those
                                     // provided on the command-line

  param.print(std::cout); // Print parameters
  std::cout.precision(15);
  const auto config = make_simulation_config(param);
  return SimulationRunner(config).run();
}
