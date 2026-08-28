#include <folded_xxz/initial_state.h>

itensor::MPS create_initial_state(const itensor::SiteSet &sites) {
  const int site_count = itensor::length(sites);
  itensor::MPS psi(sites);
  for (int first_site = 1; first_site < site_count; first_site += 2) {
    const auto left_site = sites(first_site);
    const auto right_site = sites(first_site + 1);
    itensor::ITensor wavefunction;

    if (first_site == 1) {
      const auto right_link =
          itensor::commonIndex(psi(first_site + 1), psi(first_site + 2), "Link");
      wavefunction = itensor::ITensor(left_site, right_site, right_link);
      psi.ref(first_site) = itensor::ITensor(left_site);
      psi.ref(first_site + 1) = itensor::ITensor(right_site, right_link);
      wavefunction.set(left_site(1), right_site(2), right_link(1), itensor::ISqrt2);
      wavefunction.set(left_site(2), right_site(1), right_link(1), -itensor::ISqrt2);
    } else if (first_site == site_count - 1) {
      const auto left_link = itensor::commonIndex(psi(first_site - 1), psi(first_site), "Link");
      wavefunction = itensor::ITensor(left_link, left_site, right_site);
      psi.ref(first_site) = itensor::ITensor(left_site, left_link);
      psi.ref(first_site + 1) = itensor::ITensor(right_site);
      wavefunction.set(left_link(1), left_site(1), right_site(2), itensor::ISqrt2);
      wavefunction.set(left_link(1), left_site(2), right_site(1), -itensor::ISqrt2);
    } else {
      const auto left_link = itensor::commonIndex(psi(first_site - 1), psi(first_site), "Link");
      const auto right_link =
          itensor::commonIndex(psi(first_site + 1), psi(first_site + 2), "Link");
      wavefunction = itensor::ITensor(left_link, left_site, right_site, right_link);
      psi.ref(first_site) = itensor::ITensor(left_site, left_link);
      psi.ref(first_site + 1) = itensor::ITensor(right_site, right_link);
      wavefunction.set(left_link(1), left_site(1), right_site(2), right_link(1), itensor::ISqrt2);
      wavefunction.set(left_link(1), left_site(2), right_site(1), right_link(1), -itensor::ISqrt2);
    }

    itensor::ITensor singular_values;
    itensor::svd(wavefunction, psi.ref(first_site), singular_values, psi.ref(first_site + 1));
    psi.ref(first_site) *= singular_values;
  }
  return psi;
}
