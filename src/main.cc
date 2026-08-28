#include "simulation_runner.h"

#include <folded_xxz/simulation_config.h>

#include <exception>
#include <iostream>

int main(int argc, char *argv[]) {
  try {
    return SimulationRunner(parse_simulation_config(argc, argv)).run();
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << std::endl;
    return 1;
  }
}
