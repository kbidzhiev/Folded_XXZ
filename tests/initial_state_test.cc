#include "itensor/all.h"

#include <folded_xxz/initial_state.h>

#include <cmath>
#include <cstdlib>
#include <iostream>

int main() {
  itensor::SpinHalf sites(8, {"ConserveQNs=", false});
  const auto psi = create_initial_state(sites);
  if (std::abs(itensor::innerC(psi, psi).real() - 1.0) > 1e-12) {
    std::cerr << "FAIL: initial state is not normalized" << std::endl;
    return EXIT_FAILURE;
  }
}
