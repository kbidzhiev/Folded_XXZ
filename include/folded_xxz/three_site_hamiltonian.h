#pragma once

#include "itensor/all.h"

#include <folded_xxz/model_config.h>

class ThreeSiteHamiltonian {
public:
  ThreeSiteHamiltonian(const itensor::SiteSet &sites, const ModelConfig &config);

  int center() const;
  itensor::MPO mpo() const;

private:
  void add_interaction_terms(int left, int middle, int right, double coupling);
  void add_field_terms(const ModelConfig &config);
  void add_bulk_terms(double coupling);

  const int site_count_;
  const int center_;
  itensor::AutoMPO terms_;
};
