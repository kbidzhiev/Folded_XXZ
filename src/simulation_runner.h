#pragma once

#include <folded_xxz/simulation_config.h>

class SimulationRunner {
public:
  explicit SimulationRunner(SimulationConfig config);

  int run();

private:
  SimulationConfig config_;
};
