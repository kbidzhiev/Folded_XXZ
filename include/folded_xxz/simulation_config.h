#pragma once

#include <folded_xxz/parameters.h>
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
  int site_count;
  int max_bond_dimension;
  double real_time_cutoff;
  double thermal_cutoff;
  OutputConfig output;
  SimulationSchedule schedule;
};

SimulationConfig make_simulation_config(const ThreeSiteParam &parameters);
