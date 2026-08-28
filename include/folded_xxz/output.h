#pragma once

#include "itensor/all.h"

#include <fstream>

#include <folded_xxz/simulation_config.h>
#include <folded_xxz/simulation_schedule.h>

struct OutputFiles {
  std::ofstream entropy;
  std::ofstream singular_values;
  std::ofstream entropy_profile;
  std::ofstream magnetization;
  std::ofstream average_magnetization;
  std::ofstream energy_beta;
  std::ofstream energy_profile;
  std::ofstream q1minus_profile;
  std::ofstream q2_profile;
};

OutputFiles open_output_files(const OutputConfig &config, int center);

void write_observables(itensor::MPS &psi, const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites,
                       const itensor::MPO &hamiltonian, const OutputConfig &config,
                       OutputFiles &output_files, int step, double time,
                       const SimulationSchedule &schedule, int site_count, int center);
