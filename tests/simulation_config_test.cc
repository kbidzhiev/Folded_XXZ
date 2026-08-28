#include <folded_xxz/simulation_config.h>

#include <cstdlib>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

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

SimulationConfig parse(std::initializer_list<const char *> arguments) {
  std::vector<std::string> values{"simulation_config_test"};
  values.insert(values.end(), arguments.begin(), arguments.end());
  std::vector<char *> argv;
  argv.reserve(values.size());
  for (auto &value : values)
    argv.push_back(value.data());
  return parse_simulation_config(static_cast<int>(argv.size()), argv.data());
}

} // namespace

int main() {
  const auto config = parse({"N", "4", "Entropy", "1", "SVD_spec", "0.1", "EnergyProf", "0.2"});
  expect(config.model.site_count == 8, "site count");
  expect(config.max_bond_dimension == 4000, "maximum bond dimension");
  expect(config.output.entropy, "entropy output");
  expect(config.output.singular_value_interval == 0.1, "singular-value output interval");
  expect(config.output.energy_profile_interval == 0.2, "energy-profile output interval");
  expect(config.schedule.total_steps == 20, "default schedule");

  expect_invalid_argument([&] { parse({"N", "0"}); }, "zero site count");
  expect_invalid_argument([&] { parse({"TrotterOrder", "3"}); }, "invalid Trotter order");
  expect_invalid_argument([&] { parse({"unknown", "1"}); }, "unknown parameter");
  expect_invalid_argument([&] { parse({"N"}); }, "missing parameter value");
}
