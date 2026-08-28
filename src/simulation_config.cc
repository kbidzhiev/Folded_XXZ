#include <folded_xxz/simulation_config.h>

#include <stdexcept>

SimulationConfig make_simulation_config(const ThreeSiteParam &parameters) {
  const int half_site_count = parameters.integer_value("N");
  const int max_bond_dimension = parameters.integer_value("max_bond");
  const double real_time_cutoff = parameters.value("trunc");
  const double thermal_cutoff = parameters.value("trunc0");

  if (half_site_count <= 0)
    throw std::invalid_argument("N must be positive");
  if (max_bond_dimension <= 0)
    throw std::invalid_argument("max_bond must be positive");
  if (real_time_cutoff < 0 || thermal_cutoff < 0)
    throw std::invalid_argument("trunc and trunc0 must not be negative");

  return {2 * half_site_count,
          max_bond_dimension,
          real_time_cutoff,
          thermal_cutoff,
          {parameters.value("Entropy") != 0, parameters.value("SVD_spec"),
           parameters.value("Energy_beta") > 0, parameters.value("Eprof"), parameters.value("Sz"),
           parameters.value("EnergyProf"), parameters.value("Q2Prof")},
          {parameters.value("tau"), parameters.value("dbeta"), parameters.value("T"),
           parameters.value("TL"), parameters.value("TR")}};
}
