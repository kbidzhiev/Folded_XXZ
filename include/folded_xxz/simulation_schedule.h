#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

class SimulationSchedule {
public:
  SimulationSchedule(double tau, double dbeta, double final_time, double left_temperature,
                     double right_temperature)
      : tau(tau), dbeta(dbeta) {
    if (tau <= 0 || dbeta <= 0 || left_temperature <= 0 || right_temperature <= 0) {
      throw std::invalid_argument("tau, dbeta, TL, and TR must be positive");
    }
    if (final_time < 0) {
      throw std::invalid_argument("T must not be negative");
    }
    beta_steps_min =
        steps_for(std::min(1.0 / left_temperature, 1.0 / right_temperature), "beta", dbeta);
    beta_steps_max =
        steps_for(std::max(1.0 / left_temperature, 1.0 / right_temperature), "beta", dbeta);
    real_time_steps = steps_for(final_time, "T", tau);
    total_steps = beta_steps_max + real_time_steps;
  }

  bool output_due(int step, double interval, const char *name) const {
    return step % steps_for(interval, name, tau) == 0;
  }

  const double tau;
  const double dbeta;
  int beta_steps_min;
  int beta_steps_max;
  int real_time_steps;
  int total_steps;

private:
  static int steps_for(double duration, const char *name, double step_size) {
    const double steps = duration / step_size;
    const double rounded = std::round(steps);
    if (steps < 0 || std::abs(steps - rounded) > 1e-6) {
      throw std::invalid_argument(std::string(name) + " must be a multiple of its step size");
    }
    return static_cast<int>(rounded);
  }
};
