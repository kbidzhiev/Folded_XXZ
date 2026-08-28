#include <folded_xxz/three_site_hamiltonian.h>

#include <cmath>

ThreeSiteHamiltonian::ThreeSiteHamiltonian(const itensor::SiteSet &sites, const ModelConfig &config)
    : site_count_(itensor::length(sites)), center_(site_count_ / 2), terms_(sites) {
  add_field_terms(config);
  add_bulk_terms(config.coupling);
}

int ThreeSiteHamiltonian::center() const { return center_; }

itensor::MPO ThreeSiteHamiltonian::mpo() const { return itensor::toMPO(terms_); }

void ThreeSiteHamiltonian::add_interaction_terms(int left, int middle, int right, double coupling) {
  terms_ += coupling, "S+", left, "S-", right;
  terms_ += coupling, "S-", left, "S+", right;
  terms_ += -2 * coupling, "S+", left, "Sz", middle, "S-", right;
  terms_ += -2 * coupling, "S-", left, "Sz", middle, "S+", right;
}

void ThreeSiteHamiltonian::add_field_terms(const ModelConfig &config) {
  for (int site = 1; site <= site_count_; site += 2) {
    const double field = site <= center_ ? config.left_thermal_field : config.right_thermal_field;
    terms_ += -config.coupling * field * std::pow(-1, (site + 1) / 2), "Sz", site;
  }
}

void ThreeSiteHamiltonian::add_bulk_terms(double coupling) {
  for (int left = 1; left <= site_count_ - 4; left += 2) {
    add_interaction_terms(left, left + 2, left + 4, coupling);
  }
}
