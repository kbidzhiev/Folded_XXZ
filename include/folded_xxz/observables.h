#pragma once

#include "itensor/all.h"
#include <complex>
#include <vector>

double Entropy(itensor::MPS &psi, int i, std::vector<double> &sing_vals, double r);
int BondDim(const itensor::MPS &psi, int i);
double Sz(itensor::MPS &psi, const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites, int i);
double Energy(itensor::MPS &psi, const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites, int i);
std::complex<double> Q1(itensor::MPS &psi,
                        const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites, int i);
std::complex<double> Q2(itensor::MPS &psi,
                        const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites, int i);
