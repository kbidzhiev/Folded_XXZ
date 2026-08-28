#include <folded_xxz/simulation_schedule.h>

#include <functional>
#include <iostream>
#include <stdexcept>
#include <cstdlib>

namespace {

void expect(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << std::endl;
    std::exit(1);
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
  const SimulationSchedule schedule(0.1, 0.2, 0.5, 2.5, 5.0);
  expect(schedule.beta_steps_min == 1, "minimum beta steps");
  expect(schedule.beta_steps_max == 2, "maximum beta steps");
  expect(schedule.real_time_steps == 5, "real-time steps");
  expect(schedule.total_steps == 7, "total steps");
  expect(schedule.output_due(6, 0.3, "Q2Prof"), "Q2Prof interval is due at step 6");
  expect(!schedule.output_due(5, 0.3, "Q2Prof"), "Q2Prof interval is not due at step 5");

  expect_invalid_argument([] { SimulationSchedule(0.0, 0.1, 0.0, 1.0, 1.0); }, "zero tau");
  expect_invalid_argument([] { SimulationSchedule(0.1, 0.0, 0.0, 1.0, 1.0); }, "zero dbeta");
  expect_invalid_argument([] { SimulationSchedule(0.1, 0.1, 0.0, 0.0, 1.0); }, "zero TL");
  expect_invalid_argument([] { SimulationSchedule(0.1, 0.1, 0.0, 1.0, 0.0); }, "zero TR");
  expect_invalid_argument([] { SimulationSchedule(0.1, 0.1, 0.15, 1.0, 1.0); }, "non-integral T");
  expect_invalid_argument([&] { schedule.output_due(0, 0.25, "Q2Prof"); }, "non-integral Q2Prof");
}
