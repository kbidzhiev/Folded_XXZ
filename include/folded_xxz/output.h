#pragma once

#include "itensor/all.h"

#include <memory>

#include <folded_xxz/simulation_config.h>
#include <folded_xxz/simulation_schedule.h>

class ObservableWriter {
public:
  ObservableWriter(OutputConfig config, SimulationSchedule schedule, int center);
  ~ObservableWriter();

  ObservableWriter(const ObservableWriter &) = delete;
  ObservableWriter &operator=(const ObservableWriter &) = delete;

  void write(itensor::MPS &psi, const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites,
             const itensor::MPO &hamiltonian, int step, double time);

private:
  struct Files;

  OutputConfig config_;
  SimulationSchedule schedule_;
  int center_;
  std::unique_ptr<Files> files_;
};
