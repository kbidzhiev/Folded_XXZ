#include <folded_xxz/simulation_config.h>

#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

struct RawSimulationConfig {
  double physical_site_count = 10;
  double coupling = 1.0;
  double tau = 0.01;
  double dbeta = 0.01;
  double final_time = 0;
  double left_temperature = 100;
  double right_temperature = 5;
  double left_field = 0;
  double right_field = 0;
  double entropy = 0;
  double entropy_profile_interval = 0;
  double energy_profile_interval = 0;
  double q2_profile_interval = 0;
  double magnetization_interval = 0;
  double singular_value_interval = 0;
  double max_bond_dimension = 4000;
  double real_time_cutoff = 1e-10;
  double thermal_cutoff = 1e-10;
  double trotter_order = 2;
  double energy_beta = 1;
};

double parse_double(const char *text) {
  char *end = nullptr;
  const double value = std::strtod(text, &end);
  if (end == text || *end != '\0')
    throw std::invalid_argument(std::string(text) + " is not a valid number");
  return value;
}

int as_integer(double value, const char *name) {
  const double rounded = std::round(value);
  if (std::abs(rounded - value) > 1e-6)
    throw std::invalid_argument(std::string(name) + " must be an integer");
  if (rounded < std::numeric_limits<int>::min() || rounded > std::numeric_limits<int>::max())
    throw std::out_of_range(std::string(name) + " is outside the int range");
  return static_cast<int>(rounded);
}

void set_value(RawSimulationConfig &config, const std::string &name, double value) {
  if (name == "N")
    config.physical_site_count = value;
  else if (name == "J")
    config.coupling = value;
  else if (name == "tau")
    config.tau = value;
  else if (name == "dbeta")
    config.dbeta = value;
  else if (name == "T")
    config.final_time = value;
  else if (name == "TL")
    config.left_temperature = value;
  else if (name == "TR")
    config.right_temperature = value;
  else if (name == "hL")
    config.left_field = value;
  else if (name == "hR")
    config.right_field = value;
  else if (name == "Entropy")
    config.entropy = value;
  else if (name == "Eprof")
    config.entropy_profile_interval = value;
  else if (name == "EnergyProf")
    config.energy_profile_interval = value;
  else if (name == "Q2Prof")
    config.q2_profile_interval = value;
  else if (name == "Sz")
    config.magnetization_interval = value;
  else if (name == "SVD_spec")
    config.singular_value_interval = value;
  else if (name == "max_bond")
    config.max_bond_dimension = value;
  else if (name == "trunc")
    config.real_time_cutoff = value;
  else if (name == "trunc0")
    config.thermal_cutoff = value;
  else if (name == "TrotterOrder")
    config.trotter_order = value;
  else if (name == "Energy_beta")
    config.energy_beta = value;
  else
    throw std::invalid_argument("unknown parameter: " + name);
}

} // namespace

SimulationConfig parse_simulation_config(int argc, char *argv[]) {
  RawSimulationConfig raw;
  for (int index = 1; index < argc; index += 2) {
    if (index + 1 == argc)
      throw std::invalid_argument(std::string("missing value after ") + argv[index]);
    set_value(raw, argv[index], parse_double(argv[index + 1]));
  }

  const int physical_site_count = as_integer(raw.physical_site_count, "N");
  const int max_bond_dimension = as_integer(raw.max_bond_dimension, "max_bond");
  const int trotter_order = as_integer(raw.trotter_order, "TrotterOrder");
  if (physical_site_count <= 0)
    throw std::invalid_argument("N must be positive");
  if (max_bond_dimension <= 0)
    throw std::invalid_argument("max_bond must be positive");
  if (raw.real_time_cutoff < 0 || raw.thermal_cutoff < 0)
    throw std::invalid_argument("trunc and trunc0 must not be negative");
  if (trotter_order != 1 && trotter_order != 2)
    throw std::invalid_argument("TrotterOrder must be 1 or 2");

  return {{2 * physical_site_count, raw.coupling, raw.left_field * raw.left_temperature,
           raw.right_field * raw.right_temperature, trotter_order},
          max_bond_dimension,
          raw.real_time_cutoff,
          raw.thermal_cutoff,
          {raw.entropy != 0, raw.singular_value_interval, raw.energy_beta > 0,
           raw.entropy_profile_interval, raw.magnetization_interval, raw.energy_profile_interval,
           raw.q2_profile_interval},
          {raw.tau, raw.dbeta, raw.final_time, raw.left_temperature, raw.right_temperature}};
}
