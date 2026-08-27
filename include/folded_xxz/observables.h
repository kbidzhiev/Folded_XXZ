#pragma once

#include "itensor/all.h"
#include <vector>
#include <string>
#include <complex>

double Entropy(itensor::MPS &psi, int i, std::vector<double> &sing_vals, double r);
// Bond Dim
int BondDim(const itensor::MPS &psi, int i);

// < Sz_i >
double Sz(itensor::MPS &psi, const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites, int i);

//< Sp_i Sm_i+4 >
std::complex<double> Correlation(itensor::MPS &psi,
                                 const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites,
                                 const std::string &op_name1, const std::string &op_name2, int i,
                                 int j);

std::complex<double> SzCorrelation(itensor::MPS &psi,
                                   const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites,
                                   const std::string &op_name1, const std::string &op_name2, int i);

// EgergyKin + EnergyPot at site i (i,i+2,i+4)

double Energy(itensor::MPS &psi, const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites, int i);

// Current(i,i+4) + Current_z at site i (i,i+2,i+4)
double Q1minus(itensor::MPS &psi, const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites, int i);

std::complex<double> Q1(itensor::MPS &psi,
                        const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites, int i);
std::complex<double> KKDD(itensor::MPS &psi,
                          const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites, int i);
std::complex<double> Correlations5site(itensor::MPS &psi,
                                       const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites,
                                       const std::string &op_name1, const std::string &op_name2,
                                       const std::string &op_name3, const std::string &op_name4,
                                       const std::string &op_name5, int i);
std::complex<double> Q2(itensor::MPS &psi,
                        const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites, int i);
double Q2plus(itensor::MPS &psi, const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites, int i);
double Q2minus(itensor::MPS &psi, const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites, int i);
