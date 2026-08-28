#include <folded_xxz/output.h>

#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <folded_xxz/observables.h>

struct ObservableWriter::Files {
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

ObservableWriter::ObservableWriter(OutputConfig config, SimulationSchedule schedule, int center)
    : config_(config), schedule_(schedule), center_(center), files_(std::make_unique<Files>()) {
  constexpr auto mode = std::ofstream::out;

  if (config_.entropy) {
    files_->entropy.open("Entropy_center.dat", mode);
    files_->entropy.precision(15);
    files_->entropy << "#time \t Entropy(dot) \t BondDim(dot) \t MaxBondDim\n";
  }
  if (config_.singular_value_interval > 0) {
    files_->singular_values.open("SVD_spec.dat", mode);
    files_->singular_values.precision(15);
    files_->singular_values << "#Position=" << center_ << "\t<SVD_spectrum>\t\ttime\n";
  }
  if (config_.energy_beta) {
    files_->energy_beta.open("Energy_beta.dat", mode);
    files_->energy_beta.precision(15);
    files_->energy_beta << "beta \t Energy \n";
  }
  if (config_.entropy_profile_interval > 0) {
    files_->entropy_profile.open("Entropy_profile.dat", mode);
    files_->entropy_profile.precision(15);
    files_->entropy_profile << "#Position=i-" << center_ << std::setw(16) << "\t Entropy(i)"
                            << std::setw(16)
                            << "\t Entropy_sqrt \t Entropy_state1 \t time \t\t Bond.Dim(i)\n";
  }
  if (config_.magnetization_interval > 0) {
    files_->magnetization.open("Sz_profile.dat", mode);
    files_->magnetization.precision(15);
    files_->magnetization << "#Position=i-\t<Sz_i>\t\t(-1)^i<Sz_i>\t\t\ttime\n";

    files_->average_magnetization.open("Sz_average_profile.dat", mode);
    files_->average_magnetization.precision(15);
    files_->average_magnetization << "#Position=i-\t0.5<Sz_2+1i> + 0.5<Sz_2+1i>\t\t\ttime\n";
  }
  if (config_.energy_profile_interval > 0) {
    files_->energy_profile.open("Energy_profile.dat", mode);
    files_->energy_profile.precision(15);
    files_->energy_profile << "#Position=i-\t<Ham_i>\t" << center_ << "\t\ttime(or beta)\n";

    files_->q1minus_profile.open("Q1minus_profile.dat", mode);
    files_->q1minus_profile.precision(15);
    files_->q1minus_profile << "#Position=i-\t<Q1minus_i>\t" << center_ << "\t\ttime(or beta)\n";
  }
  if (config_.q2_profile_interval > 0) {
    files_->q2_profile.open("Q2_profile.dat", mode);
    files_->q2_profile.precision(15);
    files_->q2_profile << "#Position=i-" << center_ << std::setw(16) << "\t Entropy(i)"
                       << std::setw(16) << "\t Q2plus \t Q2minus \t time \t \n";
  }
}

ObservableWriter::~ObservableWriter() = default;

void ObservableWriter::write(itensor::MPS &psi,
                             const itensor::BasicSiteSet<itensor::SpinHalfSite> &sites,
                             const itensor::MPO &hamiltonian, int step, double time) {
  const int site_count = itensor::length(sites);
  std::vector<double> singular_values;
  if (config_.entropy) {
    const double entropy = Entropy(psi, center_, singular_values, 1);
    files_->entropy << time << "\t" << std::setw(16) << std::setfill('0') << entropy << "\t"
                    << BondDim(psi, center_) << "\t" << itensor::maxLinkDim(psi) << std::endl;

    if (config_.singular_value_interval > 0 &&
        schedule_.output_due(step, config_.singular_value_interval, "SVD_spec")) {
      files_->singular_values << "\"t=" << time << "\"" << std::endl;
      for (int index = 0; index < static_cast<int>(singular_values.size()); ++index) {
        files_->singular_values << index + 1 << "\t" << singular_values[index] << "\t\t" << time
                                << std::endl;
      }
      if (step < schedule_.total_steps)
        files_->singular_values << std::endl << std::endl;
    }
  }

  if (config_.entropy_profile_interval > 0 &&
      schedule_.output_due(step, config_.entropy_profile_interval, "Eprof")) {
    files_->entropy_profile << "\"t=" << time << "\"" << std::endl;
    for (int site = 1; site < site_count; ++site) {
      const double entropy = Entropy(psi, site, singular_values, 1);
      files_->entropy_profile << site + 0.5 - center_ << "\t" << std::setw(16) << std::setfill('0')
                              << entropy << "\t" << std::setw(4) << std::setfill('0')
                              << BondDim(psi, site) << "\t" << time << std::endl;
    }
    if (step < schedule_.total_steps)
      files_->entropy_profile << "\n\n";
  }

  if (config_.energy_beta) {
    double local_energy = 0;
    double count = 0;
    const int half_width = (site_count / 4) % 2 == 0 ? site_count / 4 : site_count / 4 + 1;
    for (int site = center_ + 1 - half_width; site < center_ + 1 + half_width; site += 2) {
      local_energy += Energy(psi, sites, site);
      ++count;
    }
    local_energy /= count;

    files_->energy_beta << time << "\t"
                        << std::real(itensor::innerC(psi, hamiltonian, psi)) /
                               std::real(itensor::innerC(psi, psi))
                        << "\t" << local_energy << "\t" << itensor::maxLinkDim(psi) << std::endl;
  }

  if (config_.energy_profile_interval > 0 && schedule_.beta_steps_max <= step &&
      schedule_.output_due(step, config_.energy_profile_interval, "EnergyProf")) {
    files_->energy_profile << "\"t=" << time << "\"" << std::endl;
    files_->q1minus_profile << "\"t=" << time << "\"" << std::endl;
    for (int site = 1; site <= site_count - 5; site += 2) {
      const std::complex<double> q1 = Q1(psi, sites, site);
      files_->energy_profile << site / 2 - center_ / 2 + 1 << "\t" << std::real(q1) << "\t" << time
                             << std::endl;
      files_->q1minus_profile << site / 2 - center_ / 2 + 1 << "\t" << std::imag(q1) << "\t" << time
                              << std::endl;
    }
    files_->energy_profile << "\n\n";
    files_->q1minus_profile << "\n\n";
  }

  if (config_.q2_profile_interval > 0 && schedule_.beta_steps_max <= step &&
      schedule_.output_due(step, config_.q2_profile_interval, "Q2Prof")) {
    files_->q2_profile << "\"t=" << time << "\"" << std::endl;
    for (int site = 1; site <= site_count - 9; site += 2) {
      const std::complex<double> q2 = Q2(psi, sites, site);
      files_->q2_profile << site / 2 - center_ / 2 + 1 << "\t" << std::real(q2) << "\t"
                         << std::imag(q2) << "\t" << time << std::endl;
    }
    files_->q2_profile << "\n\n";
  }

  if (config_.magnetization_interval > 0 && schedule_.beta_steps_max <= step &&
      schedule_.output_due(step, config_.magnetization_interval, "Sz")) {
    files_->magnetization << "\"t=" << time << "\"" << std::endl;
    double total_magnetization = 0;
    double left_magnetization = 0;
    double right_magnetization = 0;
    double center_magnetization = 0;
    double previous_magnetization = 0;
    for (int site = 1; site <= site_count; site += 2) {
      const double magnetization = Sz(psi, sites, site);
      total_magnetization += magnetization;
      if (site < center_)
        left_magnetization += magnetization;
      if (site > center_)
        right_magnetization += magnetization;
      if (site == center_)
        center_magnetization += magnetization;
      files_->magnetization << site / 2 - center_ / 2 + 1 << "\t" << magnetization << "\t"
                            << std::pow(-1, (site + 1) / 2) * magnetization << "\t" << time
                            << std::endl;

      if (((site + 1) / 2) % 2 == 1) {
        previous_magnetization = magnetization;
      } else {
        files_->average_magnetization << site / 2 - center_ / 2 + 1 << "\t"
                                      << 0.5 * (magnetization + previous_magnetization) << "\t"
                                      << time << std::endl;
      }
    }
    files_->magnetization << "\n\n";
    files_->average_magnetization << "\n\n";
    std::cout << "\n<Sz_left>=" << left_magnetization << "\t<Sz_right>=" << right_magnetization
              << "\t<Sz_DOT>=" << center_magnetization << "\t<Sz_tot>=" << total_magnetization
              << std::endl;
  }
}
