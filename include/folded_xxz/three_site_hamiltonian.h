#pragma once

#include "itensor/all.h"

#include <folded_xxz/parameters.h>

class ThreeSiteHamiltonian {
public:
  ThreeSiteHamiltonian(const itensor::SiteSet &sites, const ThreeSiteParam &param);

  int center() const;
  itensor::MPO mpo() const;

private:
  void add_interaction_terms(int left, int middle, int right, double coupling);
  void add_field_terms(const ThreeSiteParam &param);
  void add_bulk_terms(double coupling);

  const int site_count_;
  const int center_;
  itensor::AutoMPO terms_;
};
