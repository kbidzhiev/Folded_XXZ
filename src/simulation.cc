#include "itensor/all.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <cmath>
#include <vector>
#include <folded_xxz/observables.h>
#include <folded_xxz/parameters.h>
#include <folded_xxz/profile.h>
#include <folded_xxz/simulation_schedule.h>
#include <folded_xxz/three_site_hamiltonian.h>

enum class EvolutionMode { ImaginaryTime, RealTime };
enum class ThermalRegion { FullChain, RightHalf };

class TrotterEvolution {
public:
  TrotterEvolution(const itensor::SiteSet &sites, const ThreeSiteParam &param,
                   std::complex<double> tau, EvolutionMode mode,
                   ThermalRegion thermal_region = ThermalRegion::FullChain)
      : site_count_(itensor::length(sites)), mode_(mode), thermal_region_(thermal_region) {
    initialize(sites, param, tau);
  }

  void evolve(itensor::MPS &psi, const itensor::Args &args) {
    for (auto &gate : gates_) {
      const int first_site = gate.first_site;
      const auto &tensor = gate.tensor;
      swap_next_sites(psi, first_site);
      swap_next_sites(psi, first_site + 3);
      psi.position(first_site + 2);
      auto wavefunction = psi(first_site + 1) * psi(first_site + 2) * psi(first_site + 3);
      wavefunction = tensor * wavefunction;
      wavefunction /= itensor::norm(wavefunction);
      wavefunction.noPrime();

      auto [left_tensor, remaining_tensor] = itensor::factor(
          wavefunction,
          {itensor::siteIndex(psi, first_site + 1), itensor::leftLinkIndex(psi, first_site + 1)},
          args);
      const auto link = itensor::commonIndex(left_tensor, remaining_tensor);
      auto [middle_tensor, right_tensor] =
          itensor::factor(remaining_tensor, {itensor::siteIndex(psi, first_site + 2), link}, args);
      psi.set(first_site + 1, left_tensor);
      psi.set(first_site + 2, middle_tensor);
      psi.set(first_site + 3, right_tensor);
      swap_next_sites(psi, first_site + 3);
      swap_next_sites(psi, first_site);
    }
  }

private:
  struct Gate {
    int first_site;
    itensor::ITensor tensor;
  };

  void initialize(const itensor::SiteSet &sites, const ThreeSiteParam &param,
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
      // Trotter gates from arxiv.org/abs/1901.04974, Eqs. (38) and (47).
      append_gate_layer(first_site, 0.5 * tau, sites, param);
      append_gate_layer(first_site + 2, 0.5 * tau, sites, param);
      append_gate_layer(first_site + 4, tau, sites, param);
      append_gate_layer(first_site + 2, 0.5 * tau, sites, param);
      append_gate_layer(first_site, 0.5 * tau, sites, param);
      return;
    }

    throw std::invalid_argument("TrotterOrder must be 1 or 2");
  }
  void append_gate_layer(int first_site, const std::complex<double> tau,
                         const itensor::SiteSet &sites, const ThreeSiteParam &param) {
    constexpr int gate_width = 4;
    constexpr int layer_spacing = 6;
    const int center = site_count_ / 2;

    for (int site = first_site; site <= site_count_ - gate_width; site += layer_spacing) {
      if (mode_ == EvolutionMode::ImaginaryTime && thermal_region_ == ThermalRegion::RightHalf &&
          site < center) {
        continue;
      }

      auto hamiltonian = local_hamiltonian(sites, site, param.value("J"));
      if (mode_ == EvolutionMode::ImaginaryTime) {
        add_thermal_field_terms(hamiltonian, sites, site, center, param);
      }
      gates_.push_back({site, itensor::expHermitian(hamiltonian, -tau)});
    }
  }

  itensor::ITensor local_hamiltonian(const itensor::SiteSet &sites, int first_site,
                                     double coupling) const {
    constexpr int physical_site_spacing = 2;
    const int middle_site = first_site + physical_site_spacing;
    const int last_site = middle_site + physical_site_spacing;
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

  void add_thermal_field_terms(itensor::ITensor &hamiltonian, const itensor::SiteSet &sites,
                               int first_site, int center, const ThreeSiteParam &param) const {
    constexpr int physical_site_spacing = 2;
    const double field = first_site < center ? param.value("hL") * param.value("TL")
                                             : param.value("hR") * param.value("TR");
    add_field_term(hamiltonian, sites, first_site, first_site, field);
    if (first_site == site_count_ - 4) {
      add_field_term(hamiltonian, sites, first_site + physical_site_spacing, first_site, field);
      add_field_term(hamiltonian, sites, first_site + 2 * physical_site_spacing, first_site, field);
    }
  }

  void add_field_term(itensor::ITensor &hamiltonian, const itensor::SiteSet &sites, int field_site,
                      int first_site, double field) const {
    constexpr int physical_site_spacing = 2;
    const int middle_site = first_site + physical_site_spacing;
    const int last_site = middle_site + physical_site_spacing;
    const double sign = std::pow(-1, (field_site + 1) / 2);
    hamiltonian += -field * sign *
                   itensor::op(sites, field_site == first_site ? "Sz" : "Id", first_site) *
                   itensor::op(sites, field_site == middle_site ? "Sz" : "Id", middle_site) *
                   itensor::op(sites, field_site == last_site ? "Sz" : "Id", last_site);
  }

  void swap_next_sites(itensor::MPS &psi, int first_site) const {
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

  const int site_count_;
  const EvolutionMode mode_;
  const ThermalRegion thermal_region_;
  std::vector<Gate> gates_;
};

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

void evolve_step(itensor::MPS &psi, int step, const SimulationSchedule &schedule, int center,
                 const itensor::Args &real_time_args, const itensor::Args &thermal_args,
                 TrotterEvolution &thermal_evolution, TrotterEvolution &half_thermal_evolution,
                 TrotterEvolution &real_time_evolution) {
  if (step < schedule.beta_steps_min) {
    std::cout << "Temperature evolution with H" << std::endl;
    thermal_evolution.evolve(psi, thermal_args);
    psi.orthogonalize(real_time_args);
    std::cout << "dot = " << center + 1 << std::endl;
    psi.normalize();
    return;
  }

  if (step < schedule.beta_steps_max) {
    std::cout << "Temperature evolution with half-chain H" << std::endl;
    std::cout << "n/beta_steps_max = " << step << "/" << schedule.beta_steps_max << std::endl;
    half_thermal_evolution.evolve(psi, thermal_args);
    psi.orthogonalize(real_time_args);
    psi.normalize();
    return;
  }

  std::cout << "Time evolution" << std::endl;
  real_time_evolution.evolve(psi, real_time_args);
  psi.orthogonalize(real_time_args);
}

struct OutputFiles {
  std::ofstream entropy;
  std::ofstream singular_values;
  std::ofstream entropy_profile;
  std::ofstream magnetization;
  std::ofstream average_magnetization;
  std::ofstream energy_beta;
  std::ofstream energy_profile;
  std::ofstream q1minus_profile;
  std::ofstream q2_profile;
};

OutputFiles open_output_files(const ThreeSiteParam &param, int center) {
  OutputFiles files;
  constexpr auto mode = std::ofstream::out;

  if (param.value("Entropy") != 0) {
    files.entropy.open("Entropy_center.dat", mode);
    files.entropy.precision(15);
    files.entropy << "#time \t Entropy(dot) \t BondDim(dot) \t MaxBondDim\n";
  }
  if (param.value("SVD_spec") > 0) {
    files.singular_values.open("SVD_spec.dat", mode);
    files.singular_values.precision(15);
    files.singular_values << "#Position=" << center << "\t<SVD_spectrum>\t\ttime\n";
  }
  if (param.value("Energy_beta") > 0) {
    files.energy_beta.open("Energy_beta.dat", mode);
    files.energy_beta.precision(15);
    files.energy_beta << "beta \t Energy \n";
  }
  if (param.value("Eprof") > 0) {
    files.entropy_profile.open("Entropy_profile.dat", mode);
    files.entropy_profile.precision(15);
    files.entropy_profile << "#Position=i-" << center << std::setw(16) << "\t Entropy(i)"
                          << std::setw(16)
                          << "\t Entropy_sqrt \t Entropy_state1 \t time \t\t Bond.Dim(i)\n";
  }
  if (param.value("Sz") > 0) {
    files.magnetization.open("Sz_profile.dat", mode);
    files.magnetization.precision(15);
    files.magnetization << "#Position=i-\t<Sz_i>\t\t(-1)^i<Sz_i>\t\t\ttime\n";

    files.average_magnetization.open("Sz_average_profile.dat", mode);
    files.average_magnetization.precision(15);
    files.average_magnetization << "#Position=i-\t0.5<Sz_2+1i> + 0.5<Sz_2+1i>\t\t\ttime\n";
  }
  if (param.value("EnergyProf") > 0) {
    files.energy_profile.open("Energy_profile.dat", mode);
    files.energy_profile.precision(15);
    files.energy_profile << "#Position=i-\t<Ham_i>\t" << center << "\t\ttime(or beta)\n";

    files.q1minus_profile.open("Q1minus_profile.dat", mode);
    files.q1minus_profile.precision(15);
    files.q1minus_profile << "#Position=i-\t<Q1minus_i>\t" << center << "\t\ttime(or beta)\n";
  }
  if (param.value("Q2Prof") > 0) {
    files.q2_profile.open("Q2_profile.dat", mode);
    files.q2_profile.precision(15);
    files.q2_profile << "#Position=i-" << center << std::setw(16) << "\t Entropy(i)"
                     << std::setw(16) << "\t Q2plus \t Q2minus \t time \t \n";
  }
  return files;
}

void write_observables(itensor::MPS &psi, const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites,
                       const itensor::MPO &hamiltonian, const ThreeSiteParam &param,
                       OutputFiles &output_files, int step, double time,
                       const SimulationSchedule &schedule, int site_count, int center) {
  std::vector<double> singular_values;
  if (param.value("Entropy") != 0) {
    const double entropy = Entropy(psi, center, singular_values, 1);
    output_files.entropy << time << "\t" << std::setw(16) << std::setfill('0') << entropy << "\t"
                         << BondDim(psi, center) << "\t" << itensor::maxLinkDim(psi) << std::endl;

    if (param.value("SVD_spec") > 0 &&
        schedule.output_due(step, param.value("SVD_spec"), "SVD_spec")) {
      output_files.singular_values << "\"t=" << time << "\"" << std::endl;
      for (int index = 0; index < static_cast<int>(singular_values.size()); ++index) {
        output_files.singular_values << index + 1 << "\t" << singular_values[index] << "\t\t"
                                     << time << std::endl;
      }
      if (step < schedule.total_steps)
        output_files.singular_values << std::endl << std::endl;
    }
  }

  if (param.value("Eprof") > 0 && schedule.output_due(step, param.value("Eprof"), "Eprof")) {
    output_files.entropy_profile << "\"t=" << time << "\"" << std::endl;
    for (int site = 1; site < site_count; ++site) {
      const double entropy = Entropy(psi, site, singular_values, 1);
      output_files.entropy_profile
          << site + 0.5 - center << "\t" << std::setw(16) << std::setfill('0') << entropy << "\t"
          << std::setw(4) << std::setfill('0') << BondDim(psi, site) << "\t" << time << std::endl;
    }
    if (step < schedule.total_steps)
      output_files.entropy_profile << "\n\n";
  }

  if (param.value("Energy_beta") > 0) {
    double local_energy = 0;
    double count = 0;
    const int half_width = (site_count / 4) % 2 == 0 ? site_count / 4 : site_count / 4 + 1;
    for (int site = center + 1 - half_width; site < center + 1 + half_width; site += 2) {
      local_energy += Energy(psi, sites, site);
      ++count;
    }
    local_energy /= count;

    output_files.energy_beta << time << "\t"
                             << std::real(itensor::innerC(psi, hamiltonian, psi)) /
                                    std::real(itensor::innerC(psi, psi))
                             << "\t" << local_energy << "\t" << itensor::maxLinkDim(psi)
                             << std::endl;
  }

  if (param.value("EnergyProf") > 0 && schedule.beta_steps_max <= step &&
      schedule.output_due(step, param.value("EnergyProf"), "EnergyProf")) {
    output_files.energy_profile << "\"t=" << time << "\"" << std::endl;
    output_files.q1minus_profile << "\"t=" << time << "\"" << std::endl;
    for (int site = 1; site <= site_count - 5; site += 2) {
      const std::complex<double> q1 = Q1(psi, sites, site);
      output_files.energy_profile << site / 2 - center / 2 + 1 << "\t" << std::real(q1) << "\t"
                                  << time << std::endl;
      output_files.q1minus_profile << site / 2 - center / 2 + 1 << "\t" << std::imag(q1) << "\t"
                                   << time << std::endl;
    }
    output_files.energy_profile << "\n\n";
    output_files.q1minus_profile << "\n\n";
  }

  if (param.value("Q2Prof") > 0 && schedule.beta_steps_max <= step &&
      schedule.output_due(step, param.value("Q2Prof"), "Q2Prof")) {
    output_files.q2_profile << "\"t=" << time << "\"" << std::endl;
    for (int site = 1; site <= site_count - 9; site += 2) {
      const std::complex<double> q2 = Q2(psi, sites, site);
      output_files.q2_profile << site / 2 - center / 2 + 1 << "\t" << std::real(q2) << "\t"
                              << std::imag(q2) << "\t" << time << std::endl;
    }
    output_files.q2_profile << "\n\n";
  }

  if (param.value("Sz") > 0 && schedule.beta_steps_max <= step &&
      schedule.output_due(step, param.value("Sz"), "Sz")) {
    output_files.magnetization << "\"t=" << time << "\"" << std::endl;
    double total_magnetization = 0;
    double left_magnetization = 0;
    double right_magnetization = 0;
    double center_magnetization = 0;
    double previous_magnetization = 0;
    for (int site = 1; site <= site_count; site += 2) {
      const double magnetization = Sz(psi, sites, site);
      total_magnetization += magnetization;
      if (site < center)
        left_magnetization += magnetization;
      if (site > center)
        right_magnetization += magnetization;
      if (site == center)
        center_magnetization += magnetization;
      output_files.magnetization << site / 2 - center / 2 + 1 << "\t" << magnetization << "\t"
                                 << std::pow(-1, (site + 1) / 2) * magnetization << "\t" << time
                                 << std::endl;

      if (((site + 1) / 2) % 2 == 1) {
        previous_magnetization = magnetization;
      } else {
        output_files.average_magnetization << site / 2 - center / 2 + 1 << "\t"
                                           << 0.5 * (magnetization + previous_magnetization) << "\t"
                                           << time << std::endl;
      }
    }
    output_files.magnetization << "\n\n";
    output_files.average_magnetization << "\n\n";
    std::cout << "\n<Sz_left>=" << left_magnetization << "\t<Sz_right>=" << right_magnetization
              << "\t<Sz_DOT>=" << center_magnetization << "\t<Sz_tot>=" << total_magnetization
              << std::endl;
  }
}

int run_simulation(int argc, char *argv[]) {
  LOG_DURATION("MAIN");
  ThreeSiteParam param;
  param.parse_arguments(argc, argv); // Now param contains the parameters, default values or those
                                     // provided on the command-line

  param.print(std::cout); // Print parameters
  std::cout.precision(15);
  const int N = 2 * param.integer_value("N");

  itensor::SpinHalf sites(N, {"ConserveQNs=", false}); // HILBERT_SPACE = itensor::SpinHalf
  std::cout << "Constructing Hamiltonian" << std::endl;
  ThreeSiteHamiltonian hamiltonian(sites, param);
  const int dot = hamiltonian.center();
  auto H = hamiltonian.mpo();
  std::cout << "Finished constructing Hamiltonian" << std::endl;
  auto psi = create_initial_state(sites);
  std::cout << "N= " << N << std::endl;
  const auto energy = itensor::inner(psi, H, psi);
  std::cout << "2. Initial energy=" << energy << " .Norm = " << itensor::inner(psi, psi)
            << std::endl;

  auto args = itensor::Args("Method=", "DensityMatrix", "Cutoff", param.value("trunc"), "MaxDim",
                            param.integer_value("max_bond"), "Normalize",
                            false); // arguments for time dynamics

  auto args0 = itensor::Args("Method=", "DensityMatrix", "Cutoff", param.value("trunc0"), "MaxDim",
                             param.integer_value("max_bond"), "Normalize",
                             false); // arguments for IMAGINARY time == initial temperature state

  auto output_files = open_output_files(param, dot);
  std::cout << "Trotter Gates for beta " << std::endl;
  const SimulationSchedule schedule(param.value("tau"), param.value("dbeta"), param.value("T"),
                                    param.value("TL"), param.value("TR"));
  TrotterEvolution expH_beta(sites, param, 0.5 * schedule.dbeta, EvolutionMode::ImaginaryTime);

  std::cout << "Trotter Gates Half for beta " << std::endl;
  TrotterEvolution expH_beta_half(sites, param, 0.5 * schedule.dbeta, EvolutionMode::ImaginaryTime,
                                  ThermalRegion::RightHalf);

  std::cout << "Trotter Gates for tau" << std::endl;
  param.set("hL", 0);
  param.set("hR", 0);
  TrotterEvolution expH(sites, param, itensor::Cplx_i * schedule.tau, EvolutionMode::RealTime);

  param.set("hL", 0);
  param.set("hR", 0);
  std::cout << "Finished preparing half-chain gates" << std::endl;

  for (int n = 0; n <= schedule.total_steps; ++n) {
    double time = (n - schedule.beta_steps_max) * schedule.tau;

    if (n < schedule.beta_steps_max) {
      time = n * schedule.dbeta;
      std::cout << "Beta(1/T) step #" << n << "/" << schedule.total_steps << "\t beta=" << time
                << std::endl;
    } else {
      std::cout << "Time step #" << n << "/" << schedule.total_steps << "\ttime=" << time
                << std::endl;
    }
    std::cout.flush();
    write_observables(psi, sites, H, param, output_files, n, time, schedule, N, dot);
    evolve_step(psi, n, schedule, dot, args, args0, expH_beta, expH_beta_half, expH);
    std::cout << "max bond dim = " << itensor::maxLinkDim(psi) << std::endl;
    std::cout << "Norm = " << std::real(itensor::innerC(psi, psi)) << std::endl;
    std::cout << "Energy = " << std::real(itensor::innerC(psi, H, psi)) << std::endl << std::endl;
  }

  std::cout << "Final observables: \n"
            << "max bond dim = " << itensor::maxLinkDim(psi) << std::endl
            << "Energy = " << std::real(itensor::innerC(psi, H, psi)) << std::endl
            << std::endl;

  std::cout << "\nTime evolution complete.\n Done ! \n";

  return 0;
}
