#include <folded_xxz/observables.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

constexpr double kTolerance = 1e-12;

void expect_near(double actual, double expected, const std::string &label) {
  if (std::abs(actual - expected) <= kTolerance)
    return;

  std::cerr << label << ": expected " << expected << ", got " << actual << '\n';
  std::exit(EXIT_FAILURE);
}

void expect_near(const std::complex<double> &actual, const std::complex<double> &expected,
                 const std::string &label) {
  expect_near(std::real(actual), std::real(expected), label + " (real)");
  expect_near(std::imag(actual), std::imag(expected), label + " (imag)");
}

itensor::MPS all_up_state(const itensor::SiteSet &sites) {
  auto state = itensor::InitState(sites);
  for (int site = 1; site <= itensor::length(sites); ++site) {
    state.set(site, "Up");
  }
  return itensor::MPS(state);
}

} // namespace

int main() {
  constexpr int kSiteCount = 10;
  auto sites = itensor::SpinHalf(kSiteCount, {"ConserveQNs=", false});
  auto psi = all_up_state(sites);

  for (int site = 1; site <= kSiteCount; ++site) {
    expect_near(Sz(psi, sites, site), 0.5, "Sz at site " + std::to_string(site));
  }

  std::vector<double> singular_values;
  expect_near(Entropy(psi, 5, singular_values, 1.0), 0.0, "product-state entropy");
  expect_near(Energy(psi, sites, 1), 0.0, "product-state energy");
  expect_near(Q1(psi, sites, 1), {0.0, 0.0}, "product-state Q1");
  expect_near(Q2(psi, sites, 1), {0.0, 0.0}, "product-state Q2");

  return EXIT_SUCCESS;
}
