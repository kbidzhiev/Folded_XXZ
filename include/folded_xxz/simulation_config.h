#pragma once

#include <folded_xxz/model_config.h>
#include <folded_xxz/simulation_schedule.h>

struct OutputConfig {
  bool entropy;
  double singular_value_interval;
  bool energy_beta;
  double entropy_profile_interval;
  double magnetization_interval;
  double energy_profile_interval;
  double q2_profile_interval;
};

struct SimulationConfig {
  ModelConfig model;
  int max_bond_dimension;
  double real_time_cutoff;
  double thermal_cutoff;
  OutputConfig output;
  SimulationSchedule schedule;
};

SimulationConfig parse_simulation_config(int argc, char *argv[]);
