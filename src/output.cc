#include <folded_xxz/output.h>

#include <complex>
#include <iomanip>
#include <iostream>

#include <folded_xxz/observables.h>

OutputFiles open_output_files(const OutputConfig &config, int center) {
  OutputFiles files;
  constexpr auto mode = std::ofstream::out;

  if (config.entropy) {
    files.entropy.open("Entropy_center.dat", mode);
    files.entropy.precision(15);
    files.entropy << "#time \t Entropy(dot) \t BondDim(dot) \t MaxBondDim\n";
  }
  if (config.singular_value_interval > 0) {
    files.singular_values.open("SVD_spec.dat", mode);
    files.singular_values.precision(15);
    files.singular_values << "#Position=" << center << "\t<SVD_spectrum>\t\ttime\n";
  }
  if (config.energy_beta) {
    files.energy_beta.open("Energy_beta.dat", mode);
    files.energy_beta.precision(15);
    files.energy_beta << "beta \t Energy \n";
  }
  if (config.entropy_profile_interval > 0) {
    files.entropy_profile.open("Entropy_profile.dat", mode);
    files.entropy_profile.precision(15);
    files.entropy_profile << "#Position=i-" << center << std::setw(16) << "\t Entropy(i)"
                          << std::setw(16)
                          << "\t Entropy_sqrt \t Entropy_state1 \t time \t\t Bond.Dim(i)\n";
  }
  if (config.magnetization_interval > 0) {
    files.magnetization.open("Sz_profile.dat", mode);
    files.magnetization.precision(15);
    files.magnetization << "#Position=i-\t<Sz_i>\t\t(-1)^i<Sz_i>\t\t\ttime\n";

    files.average_magnetization.open("Sz_average_profile.dat", mode);
    files.average_magnetization.precision(15);
    files.average_magnetization << "#Position=i-\t0.5<Sz_2+1i> + 0.5<Sz_2+1i>\t\t\ttime\n";
  }
  if (config.energy_profile_interval > 0) {
    files.energy_profile.open("Energy_profile.dat", mode);
    files.energy_profile.precision(15);
    files.energy_profile << "#Position=i-\t<Ham_i>\t" << center << "\t\ttime(or beta)\n";

    files.q1minus_profile.open("Q1minus_profile.dat", mode);
    files.q1minus_profile.precision(15);
    files.q1minus_profile << "#Position=i-\t<Q1minus_i>\t" << center << "\t\ttime(or beta)\n";
  }
  if (config.q2_profile_interval > 0) {
    files.q2_profile.open("Q2_profile.dat", mode);
    files.q2_profile.precision(15);
    files.q2_profile << "#Position=i-" << center << std::setw(16) << "\t Entropy(i)"
                     << std::setw(16) << "\t Q2plus \t Q2minus \t time \t \n";
  }
  return files;
}

void write_observables(itensor::MPS &psi, const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites,
                       const itensor::MPO &hamiltonian, const OutputConfig &config,
                       OutputFiles &output_files, int step, double time,
                       const SimulationSchedule &schedule, int site_count, int center) {
  std::vector<double> singular_values;
  if (config.entropy) {
    const double entropy = Entropy(psi, center, singular_values, 1);
    output_files.entropy << time << "\t" << std::setw(16) << std::setfill('0') << entropy << "\t"
                         << BondDim(psi, center) << "\t" << itensor::maxLinkDim(psi) << std::endl;

    if (config.singular_value_interval > 0 &&
        schedule.output_due(step, config.singular_value_interval, "SVD_spec")) {
      output_files.singular_values << "\"t=" << time << "\"" << std::endl;
      for (int index = 0; index < static_cast<int>(singular_values.size()); ++index) {
        output_files.singular_values << index + 1 << "\t" << singular_values[index] << "\t\t"
                                     << time << std::endl;
      }
      if (step < schedule.total_steps)
        output_files.singular_values << std::endl << std::endl;
    }
  }

  if (config.entropy_profile_interval > 0 &&
      schedule.output_due(step, config.entropy_profile_interval, "Eprof")) {
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

  if (config.energy_beta) {
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

  if (config.energy_profile_interval > 0 && schedule.beta_steps_max <= step &&
      schedule.output_due(step, config.energy_profile_interval, "EnergyProf")) {
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

  if (config.q2_profile_interval > 0 && schedule.beta_steps_max <= step &&
      schedule.output_due(step, config.q2_profile_interval, "Q2Prof")) {
    output_files.q2_profile << "\"t=" << time << "\"" << std::endl;
    for (int site = 1; site <= site_count - 9; site += 2) {
      const std::complex<double> q2 = Q2(psi, sites, site);
      output_files.q2_profile << site / 2 - center / 2 + 1 << "\t" << std::real(q2) << "\t"
                              << std::imag(q2) << "\t" << time << std::endl;
    }
    output_files.q2_profile << "\n\n";
  }

  if (config.magnetization_interval > 0 && schedule.beta_steps_max <= step &&
      schedule.output_due(step, config.magnetization_interval, "Sz")) {
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
