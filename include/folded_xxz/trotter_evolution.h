#pragma once

#include "itensor/all.h"

#include <complex>
#include <vector>

#include <folded_xxz/parameters.h>

enum class EvolutionMode { ImaginaryTime, RealTime };
enum class ThermalRegion { FullChain, RightHalf };

class TrotterEvolution {
public:
  TrotterEvolution(const itensor::SiteSet &sites, const ThreeSiteParam &param,
                   std::complex<double> tau, EvolutionMode mode,
                   ThermalRegion thermal_region = ThermalRegion::FullChain);

  void evolve(itensor::MPS &psi, const itensor::Args &args) const;
  std::size_t gate_count() const;

private:
  struct Gate {
    int first_site;
    itensor::ITensor tensor;
  };

  void build_gates(const itensor::SiteSet &sites, const ThreeSiteParam &param,
                   std::complex<double> tau);
  void apply_gate(itensor::MPS &psi, const Gate &gate, const itensor::Args &args) const;
  void append_gate_layer(int first_site, std::complex<double> tau, const itensor::SiteSet &sites,
                         const ThreeSiteParam &param);
  itensor::ITensor local_hamiltonian(const itensor::SiteSet &sites, int first_site,
                                     double coupling) const;
  void add_thermal_field_terms(itensor::ITensor &hamiltonian, const itensor::SiteSet &sites,
                               int first_site, int center, const ThreeSiteParam &param) const;
  void add_field_term(itensor::ITensor &hamiltonian, const itensor::SiteSet &sites, int field_site,
                      int first_site, double field) const;
  void swap_next_sites(itensor::MPS &psi, int first_site) const;

  const int site_count_;
  const EvolutionMode mode_;
  const ThermalRegion thermal_region_;
  std::vector<Gate> gates_;
};
