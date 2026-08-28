#include "itensor/all.h"
#include <folded_xxz/observables.h>

double Entropy(itensor::MPS &psi, const int i, std::vector<double> &sing_vals, const double r) {
  const auto bond_index = itensor::rightLinkIndex(psi, i);
  const int bond_dim = itensor::dim(bond_index);
  psi.position(i);
  itensor::ITensor wf = psi(i) * psi(i + 1);
  auto U = psi(i);
  itensor::ITensor S, V;
  const auto spectrum = itensor::svd(wf, U, S, V, {"MaxDim", bond_dim});
  itensor::Real SvN = 0.;
  sing_vals.resize(spectrum.numEigsKept());
  int j = 0;
  for (double p : spectrum.eigs()) {
    if (p > 1e-15) {
      sing_vals[j] = p;
      j++;
      SvN += -std::pow(p, r) * r * std::log(p);
    }
  }
  return SvN;
}

int BondDim(const itensor::MPS &psi, const int i) {
  const auto bond_index = itensor::rightLinkIndex(psi, i);
  return itensor::dim(bond_index);
}

double Sz(itensor::MPS &psi, const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites,
          const int i) {
  psi.position(i);
  auto ket = psi(i);
  ket *= itensor::op(sites, "Sz", i);
  ket *= itensor::dag(itensor::prime(psi(i), "Site"));
  return std::real(itensor::eltC(ket));
}

namespace {

std::complex<double> Correlation(itensor::MPS &psi,
                                 const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites,
                                 const std::string &op_name1, const std::string &op_name2,
                                 const int i, const int j) {
  itensor::ITensor ket = psi(i);
  auto Sp = itensor::op(sites, op_name1, i);
  auto Sm = itensor::op(sites, op_name2, j);
  ket *= Sp;
  auto ir = itensor::commonIndex(psi(i), psi(i + 1), "Link");
  ket *= itensor::dag(itensor::prime(itensor::prime(psi(i), "Site"), ir));
  for (int k = i + 1; k < j; ++k) {
    ket *= psi(k);
    auto right = itensor::commonIndex(psi(k), psi(k + 1), "Link");
    auto left = itensor::commonIndex(psi(k - 1), psi(k), "Link");
    ket *= itensor::dag(itensor::prime(itensor::prime(psi(k), right), left));
  }
  ket *= psi(j);
  ket *= Sm;
  auto il = itensor::commonIndex(psi(j), psi(j - 1), "Link");
  ket *= itensor::dag(itensor::prime(itensor::prime(psi(j), "Site"), il));
  std::complex<double> correlation = itensor::eltC(ket);
  return correlation;
}

std::complex<double> SzCorrelation(itensor::MPS &psi,
                                   const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites,
                                   const std::string &op_name1, const std::string &op_name2,
                                   const int i) {
  itensor::ITensor ket = psi(i);
  auto Sp = itensor::op(sites, op_name1, i);
  auto Sm = itensor::op(sites, op_name2, i + 4);
  auto Sz = itensor::op(sites, "Sz", i + 2);
  ket *= Sp;
  const auto ir = itensor::commonIndex(psi(i), psi(i + 1), "Link");
  ket *= itensor::dag(itensor::prime(itensor::prime(psi(i), "Site"), ir));
  ket *= psi(i + 1);
  ket *= itensor::dag(itensor::prime(psi(i + 1), "Link"));
  ket *= psi(i + 2);
  ket *= Sz;
  ket *= itensor::dag(itensor::prime(itensor::prime(psi(i + 2), "Site"), "Link"));
  ket *= psi(i + 3);
  ket *= itensor::dag(itensor::prime(psi(i + 3), "Link"));
  ket *= psi(i + 4);
  ket *= Sm;
  auto il = itensor::commonIndex(psi(i + 4), psi(i + 3), "Link");
  ket *= itensor::dag(itensor::prime(itensor::prime(psi(i + 4), "Site"), il));
  std::complex<double> correlation = itensor::eltC(ket);
  return correlation;
}

} // namespace

double Energy(itensor::MPS &psi, const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites,
              const int i) {
  psi.position(i);
  /*
     the "strange" coefficients 4 and 8 here appeared cause we use
     sigma matrices in the paper and spin-1/2 in the code.
     sigma_j = 2 * spin_matrix_j, so
     4 comes from XX+YY term
     8 comes from (XX+YY)*Z

     coeff 0.25 appeared cause the Hamiltonian H = 1/2 * (XX+YY)(1-Z)
     then (XX+YY) = 0.5(SpSm+SmSp), so
     H = 0.25(SpSm+SmSp)(1-Sz)

   */
  double energy_kin = 2 * std::real(4 * 0.5 * Correlation(psi, sites, "S+", "S-", i, i + 4));
  double energy_pot = 2 * std::real(-8 * 0.5 * SzCorrelation(psi, sites, "S+", "S-", i));
  double energy = 0.5 * (energy_kin + energy_pot);
  return energy;
}

std::complex<double> Q1(itensor::MPS &psi,
                        const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites, const int i) {
  psi.position(i);
  std::complex<double> q_kin = 2 * 4 * 0.5 * Correlation(psi, sites, "S-", "S+", i, i + 4);
  std::complex<double> q_pot = 2 * (-8) * 0.5 * SzCorrelation(psi, sites, "S-", "S+", i);
  const double q1_plus = 0.5 * std::real(q_kin + q_pot);
  double q1_minus = -0.5 * std::imag(q_kin + q_pot);

  return {q1_plus, q1_minus};
}

namespace {

std::complex<double> KKDD(itensor::MPS &psi,
                          const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites, const int i) {
  psi.position(i + 2);
  itensor::ITensor ket = psi(i + 2);
  auto Sp1 = itensor::op(sites, "S-", i + 2);
  auto Sp2 = itensor::op(sites, "S-", i + 4);
  auto Sm1 = itensor::op(sites, "S+", i + 6);
  auto Sm2 = itensor::op(sites, "S+", i + 8);
  ket *= Sp1;
  const auto ir = itensor::commonIndex(psi(i + 2), psi(i + 3), "Link");
  ket *= itensor::dag(itensor::prime(itensor::prime(psi(i + 2), "Site"), ir));
  ket *= psi(i + 3);
  ket *= itensor::dag(itensor::prime(psi(i + 3), "Link"));
  ket *= psi(i + 4);
  ket *= Sp2;
  ket *= itensor::dag(itensor::prime(itensor::prime(psi(i + 4), "Site"), "Link"));
  ket *= psi(i + 5);
  ket *= itensor::dag(itensor::prime(psi(i + 5), "Link"));
  ket *= psi(i + 6);
  ket *= Sm1;
  ket *= itensor::dag(itensor::prime(itensor::prime(psi(i + 6), "Site"), "Link"));
  ket *= psi(i + 7);
  ket *= itensor::dag(itensor::prime(psi(i + 7), "Link"));
  ket *= psi(i + 8);
  ket *= Sm2;
  auto il = itensor::commonIndex(psi(i + 7), psi(i + 8), "Link");
  ket *= itensor::dag(itensor::prime(itensor::prime(psi(i + 8), "Site"), il));
  const std::complex<double> kkdd = 0.25 * 16 * 2 * itensor::eltC(ket);
  return kkdd;
}

std::complex<double> Correlations5site(itensor::MPS &psi,
                                       const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites,
                                       const std::string &op_name1, const std::string &op_name2,
                                       const std::string &op_name3, const std::string &op_name4,
                                       const std::string &op_name5, const int i) {
  auto Op0 = itensor::op(sites, op_name1, i);
  auto Op2 = itensor::op(sites, op_name2, i + 2);
  auto Op4 = itensor::op(sites, op_name3, i + 4);
  auto Op6 = itensor::op(sites, op_name4, i + 6);
  auto Op8 = itensor::op(sites, op_name5, i + 8);

  itensor::ITensor ket = psi(i);
  ket *= Op0;
  const auto ir = itensor::commonIndex(psi(i), psi(i + 1), "Link");
  ket *= itensor::dag(itensor::prime(itensor::prime(psi(i), "Site"), ir));

  ket *= psi(i + 1);
  ket *= itensor::dag(itensor::prime(psi(i + 1), "Link"));

  ket *= psi(i + 2);
  ket *= Op2;
  ket *= itensor::dag(itensor::prime(itensor::prime(psi(i + 2), "Site"), "Link"));

  ket *= psi(i + 3);
  ket *= itensor::dag(itensor::prime(psi(i + 3), "Link"));

  ket *= psi(i + 4);
  ket *= Op4;
  ket *= itensor::dag(itensor::prime(itensor::prime(psi(i + 4), "Site"), "Link"));

  ket *= psi(i + 5);
  ket *= itensor::dag(itensor::prime(psi(i + 5), "Link"));

  ket *= psi(i + 6);
  ket *= Op6;
  ket *= itensor::dag(itensor::prime(itensor::prime(psi(i + 6), "Site"), "Link"));

  ket *= psi(i + 7);
  ket *= itensor::dag(itensor::prime(psi(i + 7), "Link"));

  ket *= psi(i + 8);
  ket *= Op8;
  auto il = itensor::commonIndex(psi(i + 7), psi(i + 8), "Link");
  ket *= itensor::dag(itensor::prime(itensor::prime(psi(i + 8), "Site"), il));

  std::complex<double> correlation = itensor::eltC(ket);
  return correlation;
}

} // namespace

std::complex<double> Q2(itensor::MPS &psi,
                        const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites, const int i) {
  std::complex<double> kkdd = KKDD(psi, sites, i);
  psi.position(i);
  std::complex<double> zk =
      0.5 * 8 * (Correlations5site(psi, sites, "S-", "Id", "Sz", "Id", "S+", i));
  std::complex<double> zzzk =
      0.5 * 32 * (Correlations5site(psi, sites, "S-", "Sz", "Sz", "Sz", "S+", i));
  std::complex<double> zzk1 =
      0.5 * 16 * (Correlations5site(psi, sites, "S-", "Sz", "Sz", "Id", "S+", i));
  std::complex<double> zzk2 =
      0.5 * 16 * (Correlations5site(psi, sites, "S-", "Id", "Sz", "Sz", "S+", i));
  std::complex<double> q2 = 0.25 * (kkdd + zk + zzzk - zzk1 - zzk2);
  double q2plus = -2 * std::real(q2);
  double q2minus = 2 * std::imag(q2);
  return {q2plus, q2minus};
}
