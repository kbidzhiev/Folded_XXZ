#include <folded_xxz/trotter_evolution.h>

#include <cmath>
#include <iostream>
#include <stdexcept>

TrotterEvolution::TrotterEvolution(const itensor::SiteSet &sites, const ThreeSiteParam &param,
                                   std::complex<double> tau, EvolutionMode mode,
                                   ThermalRegion thermal_region)
    : site_count_(itensor::length(sites)), mode_(mode), thermal_region_(thermal_region) {
  build_gates(sites, param, tau);
}

void TrotterEvolution::evolve(itensor::MPS &psi, const itensor::Args &args) const {
  for (const auto &gate : gates_)
    apply_gate(psi, gate, args);
}

std::size_t TrotterEvolution::gate_count() const { return gates_.size(); }

void TrotterEvolution::build_gates(const itensor::SiteSet &sites, const ThreeSiteParam &param,
                                   const std::complex<double> tau) {
  const int first_site = 1;
  const int order = param.integer_value("TrotterOrder");
  if (order == 1) {
    std::cout << "First-order Trotter scheme" << std::endl;
    append_gate_layer(first_site, tau, sites, param);
    append_gate_layer(first_site + 2, tau, sites, param);
    append_gate_layer(first_site + 4, tau, sites, param);
    return;
  }
  if (order == 2) {
    std::cout << "Second-order Trotter scheme" << std::endl;
    append_gate_layer(first_site, 0.5 * tau, sites, param);
    append_gate_layer(first_site + 2, 0.5 * tau, sites, param);
    append_gate_layer(first_site + 4, tau, sites, param);
    append_gate_layer(first_site + 2, 0.5 * tau, sites, param);
    append_gate_layer(first_site, 0.5 * tau, sites, param);
    return;
  }
  throw std::invalid_argument("TrotterOrder must be 1 or 2");
}

void TrotterEvolution::apply_gate(itensor::MPS &psi, const Gate &gate,
                                  const itensor::Args &args) const {
  const int first_site = gate.first_site;
  swap_next_sites(psi, first_site);
  swap_next_sites(psi, first_site + 3);
  psi.position(first_site + 2);
  auto wavefunction = psi(first_site + 1) * psi(first_site + 2) * psi(first_site + 3);
  wavefunction = gate.tensor * wavefunction;
  wavefunction /= itensor::norm(wavefunction);
  wavefunction.noPrime();
  auto [left_tensor, remaining_tensor] = itensor::factor(
      wavefunction,
      {itensor::siteIndex(psi, first_site + 1), itensor::leftLinkIndex(psi, first_site + 1)}, args);
  const auto link = itensor::commonIndex(left_tensor, remaining_tensor);
  auto [middle_tensor, right_tensor] =
      itensor::factor(remaining_tensor, {itensor::siteIndex(psi, first_site + 2), link}, args);
  psi.set(first_site + 1, left_tensor);
  psi.set(first_site + 2, middle_tensor);
  psi.set(first_site + 3, right_tensor);
  swap_next_sites(psi, first_site + 3);
  swap_next_sites(psi, first_site);
}

void TrotterEvolution::append_gate_layer(int first_site, const std::complex<double> tau,
                                         const itensor::SiteSet &sites,
                                         const ThreeSiteParam &param) {
  constexpr int gate_width = 4;
  constexpr int layer_spacing = 6;
  const int center = site_count_ / 2;
  for (int site = first_site; site <= site_count_ - gate_width; site += layer_spacing) {
    if (mode_ == EvolutionMode::ImaginaryTime && thermal_region_ == ThermalRegion::RightHalf &&
        site < center)
      continue;
    const double coupling = param.value("J");
    auto hamiltonian = local_hamiltonian(sites, site, coupling);
    if (mode_ == EvolutionMode::ImaginaryTime)
      add_thermal_field_terms(hamiltonian, sites, site, center, coupling, param);
    gates_.push_back({site, itensor::expHermitian(hamiltonian, -tau)});
  }
}

itensor::ITensor TrotterEvolution::local_hamiltonian(const itensor::SiteSet &sites, int first_site,
                                                     double coupling) const {
  const int middle_site = first_site + 2;
  const int last_site = middle_site + 2;
  auto hamiltonian = coupling * itensor::op(sites, "Sp", first_site) *
                     itensor::op(sites, "Id", middle_site) * itensor::op(sites, "Sm", last_site);
  hamiltonian += coupling * itensor::op(sites, "Sm", first_site) *
                 itensor::op(sites, "Id", middle_site) * itensor::op(sites, "Sp", last_site);
  hamiltonian += -2 * coupling * itensor::op(sites, "Sp", first_site) *
                 itensor::op(sites, "Sz", middle_site) * itensor::op(sites, "Sm", last_site);
  hamiltonian += -2 * coupling * itensor::op(sites, "Sm", first_site) *
                 itensor::op(sites, "Sz", middle_site) * itensor::op(sites, "Sp", last_site);
  return hamiltonian;
}

void TrotterEvolution::add_thermal_field_terms(itensor::ITensor &hamiltonian,
                                               const itensor::SiteSet &sites, int first_site,
                                               int center, double coupling,
                                               const ThreeSiteParam &param) const {
  const double field = first_site < center ? param.value("hL") * param.value("TL")
                                           : param.value("hR") * param.value("TR");
  add_field_term(hamiltonian, sites, first_site, first_site, coupling, field);
  if (first_site == site_count_ - 4) {
    add_field_term(hamiltonian, sites, first_site + 2, first_site, coupling, field);
    add_field_term(hamiltonian, sites, first_site + 4, first_site, coupling, field);
  }
}

void TrotterEvolution::add_field_term(itensor::ITensor &hamiltonian, const itensor::SiteSet &sites,
                                      int field_site, int first_site, double coupling,
                                      double field) const {
  const int middle_site = first_site + 2;
  const int last_site = middle_site + 2;
  const double sign = std::pow(-1, (field_site + 1) / 2);
  hamiltonian += -coupling * field * sign *
                 itensor::op(sites, field_site == first_site ? "Sz" : "Id", first_site) *
                 itensor::op(sites, field_site == middle_site ? "Sz" : "Id", middle_site) *
                 itensor::op(sites, field_site == last_site ? "Sz" : "Id", last_site);
}

void TrotterEvolution::swap_next_sites(itensor::MPS &psi, int first_site) const {
  psi.position(first_site);
  auto wavefunction = psi(first_site) * psi(first_site + 1);
  auto [left_tensor, right_tensor] = itensor::factor(
      wavefunction,
      {itensor::siteIndex(psi, first_site + 1), itensor::leftLinkIndex(psi, first_site)},
      {"Truncate=", false});
  psi.set(first_site, left_tensor);
  psi.set(first_site + 1, right_tensor);
  psi.position(first_site);
}
