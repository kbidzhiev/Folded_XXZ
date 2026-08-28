#include <folded_xxz/simulation_config.h>

#include <cstdlib>
#include <functional>
#include <iostream>
#include <stdexcept>

namespace {

void expect(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void expect_invalid_argument(const std::function<void()> &operation, const char *message) {
  try {
    operation();
  } catch (const std::invalid_argument &) {
    return;
  }
  expect(false, message);
}

} // namespace

int main() {
  ThreeSiteParam parameters;
  parameters.set("N", 4);
  parameters.set("Entropy", 1);
  parameters.set("SVD_spec", 0.1);
  parameters.set("EnergyProf", 0.2);

  const auto config = make_simulation_config(parameters);
  expect(config.site_count == 8, "site count");
  expect(config.max_bond_dimension == 4000, "maximum bond dimension");
  expect(config.output.entropy, "entropy output");
  expect(config.output.singular_value_interval == 0.1, "singular-value output interval");
  expect(config.output.energy_profile_interval == 0.2, "energy-profile output interval");
  expect(config.schedule.total_steps == 20, "default schedule");

  parameters.set("N", 0);
  expect_invalid_argument([&] { make_simulation_config(parameters); }, "zero site count");
}
