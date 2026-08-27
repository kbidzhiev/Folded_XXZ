#include "itensor/all.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>
#include <limits>
#include <sstream>
#include <complex>
#include <cstdlib>
#include <math.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <exception>
#include <folded_xxz/observables.h>
#include <folded_xxz/profile.h>

template <class T> std::string to_string(const T &t, unsigned int precision = 0) {
  std::stringstream ss;
  if (precision > 0)
    ss.precision(precision);
  ss << t;
  return ss.str();
}

double char2double(char *a) {
  char *end_ptr;
  const double x = std::strtod(a, &end_ptr);
  if (end_ptr == a || ('\0' != *end_ptr))
    std::cout << std::endl
              << "ERROR :" << a << " is not a valid format for a double." << std::endl,
        std::exit(0);
  return x;
}

class Parameters {
public:
  double value(const std::string &name) const {
    auto it = values_.find(name);
    if (it == values_.end())
      throw std::out_of_range("unknown parameter: " + name);
    return it->second;
  }
  int integer_value(const std::string &name) const {
    double v = value(name);
    const double rounded = std::round(v);
    if (std::abs(rounded - v) > 1e-6) {
      throw std::invalid_argument("parameter " + name + " must be an integer");
    }
    if (rounded < std::numeric_limits<int>::min() || rounded > std::numeric_limits<int>::max()) {
      throw std::out_of_range("parameter " + name + " is outside the int range");
    }
    return static_cast<int>(rounded);
  }
  void set(const std::string &name, double value) {
    if (!contains(name))
      throw std::out_of_range("unknown parameter: " + name);
    values_.at(name) = value;
  }
  bool contains(const std::string &name) const { return values_.find(name) != values_.end(); }
  void print(std::ostream &os) const {
    for (const auto &[name, value] : values_) {
      os << name << "=" << value << std::endl;
    }
  }
  void parse_arguments(int argc, char *argv[]) {
    for (int n = 1; n < argc; n++) {
      std::string var_name(argv[n]);
      if (!contains(var_name)) {
        std::cerr << "Syntax error :" << var_name << std::endl;
        std::cout << "List of command-line parameters :";
        print(std::cout);
        std::exit(0);
      }

      n++;
      if (n == argc)
        std::cerr << "Error: missing value after " << var_name << std::endl, std::exit(0);
      set(var_name, char2double(argv[n]));
    }
  }

protected:
  void define(const std::string &name, double value) { values_.emplace(name, value); }

private:
  std::map<std::string, double> values_;
};

class ThreeSiteParam : public Parameters {
public:
  ThreeSiteParam() {
    define("N", 10); // Length of the chain N-n
    define("begin", 1);
    define("end", 10);
    define("J", 1.0);
    define("tau", 0.01); // time step for the unitary evolution
    define("dbeta", 0.01);
    define("T", 0); // Total (final) time
    define("TL", 100);
    define("TR", 5);
    define("hL", 0); // alternating chemical potential or staggered magnetization
    define("hR", 0);
    define("Entropy", 0); // entanglement entropy p*log*p between left and right parts of system
    define("Eprof",
           0); // Entropy profile - parameter 0 -> nothing, dt>0 each second=integer parameter
    define("EnergyProf", 0);
    define("Q2Prof", 0);
    define("CurrentProf", 0);
    define("Current", 0);
    define("Sz", 0);
    define("SVD_spec", 0);    // SVD spectrum
    define("max_bond", 4000); // maximum bond dimension
    define("trunc", 1e-10);   // maximum truncation error
    define("trunc0", 1e-10);
    define("energy", 1e-10); // energy convergence criterion
    define("sweeps", 999);   // maximum number of sweeps in the DMRG
    define("TrotterOrder", 2);
    define("antal", 0);
    define("XXZ", 0);
    define("PBC", 0);
    define("beta", 1);
    define("Energy_beta", 1);
    define("write_wf", 0); // 0->do not write the w.-f. to disk. if dt>0 => w.-f. to disk
                           // (over)written to disk every time=dt*integer.
  }
};

class ThreeSiteHamiltonian {
public:
  int dot;
  itensor::AutoMPO ampo;
  ThreeSiteHamiltonian(const itensor::SiteSet &sites, const ThreeSiteParam &param) : ampo(sites) {
    N = itensor::length(sites);
    init(param);
    dot = N / 2;
    std::cout << "A Hamiltonian with " << N << " sites was constructed." << std::endl;
  }

private:
  int N;
  void init(const ThreeSiteParam &param) { //.init (param)
    const double J = param.value("J");
    const double n = param.value("begin");
    const double hL = param.value("hL");
    const double hR = param.value("hR");
    const double TL = param.value("TL");
    const double TR = param.value("TR");
    for (int j = n; j <= N; j += 2) {
      const double mu = j <= N / 2 ? hL * TL : hR * TR;
      ampo += -J * (mu)*std::pow(-1, (j + 1) / 2), "Sz", j;
    }
    for (int j = n; j <= N - 4; j += 2) {
      // Coefficients convert spin-1/2 operators to Pauli-matrix conventions.
      ampo += J, "S+", j, "S-", j + 4;
      ampo += J, "S-", j, "S+", j + 4;
      ampo += -2 * J, "S+", j, "Sz", j + 2, "S-", j + 4;
      ampo += -2 * J, "S-", j, "Sz", j + 2, "S+", j + 4;
    }

    std::cout << "Three-site Hamiltonian constructed" << std::endl;
    if (param.value("PBC")) {
      // This part realizes Periodic Boundary Condition (PBC)
      // term (N-1,N,1)
      ampo += J, "S+", N - 3, "S-", 1;
      ampo += J, "S-", N - 3, "S+", 1;
      ampo += -2 * J, "S+", N - 3, "Sz", N - 1, "S-", 1;
      ampo += -2 * J, "S-", N - 3, "Sz", N - 1, "S+", 1;
      std::cout << "PBC;\n sites (" << (N - 3 + 1) / 2 << ", " << (N - 1 + 1) / 2 << ", "
                << (1 + 1) / 2 << "), ";
      // term (N,1,2)
      ampo += J, "S+", N - 1, "S-", 3;
      ampo += J, "S-", N - 1, "S+", 3;
      ampo += -2 * J, "S+", N - 1, "Sz", 1, "S-", 3;
      ampo += -2 * J, "S-", N - 1, "Sz", 1, "S+", 3;
      std::cout << "(" << (N - 1 + 2) / 2 << ", " << (1 + 1) / 2 << ", " << (3 + 1) / 2 << ")"
                << std::endl;
      std::cout << "Three-site Hamiltonian uses periodic boundaries" << std::endl;
    }
  }
};

class TrotterExp {
public:
  struct TGate {
    int i1 = 0;
    itensor::ITensor G;
    TGate() {}
    TGate(int i1_, itensor::ITensor G_) : i1(i1_), G(G_) {}
  };
  TrotterExp(const itensor::SiteSet &sites, const ThreeSiteParam &param,
             const std::complex<double> tau) {
    initialize(sites, param, tau);
  };
  void initialize(const itensor::SiteSet &sites, const ThreeSiteParam &param,
                  const std::complex<double> tau) {
    const int begin = 1;
    const int end = param.value("end");
    const int order = param.value("TrotterOrder");
    if (order == 1) {
      std::cout << "First-order Trotter scheme" << std::endl;
      if (std::imag(tau) < 1e-8) {
        std::cout << "Temperature evolution" << std::endl;
        TemperatureGates(begin, end, tau, sites, param);
        TemperatureGates(begin + 2, end, tau, sites, param);
        TemperatureGates(begin + 4, end, tau, sites, param);
      } else {
        TimeGates(begin, end, tau, sites, param);
        TimeGates(begin + 2, end, tau, sites, param);
        TimeGates(begin + 4, end, tau, sites, param);
      }
    } else {
      std::cout << "Second-order Trotter scheme" << std::endl;
      double begin0 = begin;
      double begin2 = begin + 2;
      double begin4 = begin + 4;
      // Trotter gates from arxiv.org/abs/1901.04974, Eqs. (38) and (47).
      if (std::imag(tau) < 1e-8) {
        std::cout << "Temperature evolution" << std::endl;
        TemperatureGates(begin0, end, 0.5 * tau, sites, param); // A
        TemperatureGates(begin2, end, 0.5 * tau, sites, param); // B
        TemperatureGates(begin4, end, tau, sites, param);       // C
        TemperatureGates(begin2, end, 0.5 * tau, sites, param); // B
        TemperatureGates(begin0, end, 0.5 * tau, sites, param); // A
      } else {
        std::cout << "Time evolution" << std::endl;
        TimeGates(begin0, end, 0.5 * tau, sites, param); // A
        TimeGates(begin2, end, 0.5 * tau, sites, param); // B
        TimeGates(begin4, end, tau, sites, param);       // C
        TimeGates(begin2, end, 0.5 * tau, sites, param); // B
        TimeGates(begin0, end, 0.5 * tau, sites, param); // A
      }
    }
  }
  void TemperatureGates(const int begin, const int end, const std::complex<double> tau,
                        const itensor::SiteSet &sites, const ThreeSiteParam &param) {
    const int step = 6;
    const double J = param.value("J");
    const double hL = param.value("hL");
    const double hR = param.value("hR");
    const double TL = param.value("TL");
    const double TR = param.value("TR");
    const int h_half = param.value("begin");
    const int dot = itensor::length(sites) / 2;
    std::cout << "Dot in Trotter evolution = " << dot << std::endl;
    for (int j = begin; j <= end - 5; j += step) {
      if (h_half > 5 && j < dot) { //&& dot < j + step - 1
        std::cout << "j = [" << j << ", " << j + 2 << ", " << j + 4 << "]" << std::endl;
        continue; // Skip this gate and continue with the next triplet.
      }
      std::cout << "j = (" << j << ", " << j + 2 << ", " << j + 4 << ")" << std::endl;
      auto hh = J * itensor::op(sites, "Sp", j) * itensor::op(sites, "Id", j + 2) *
                itensor::op(sites, "Sm", j + 4);

      hh += J * itensor::op(sites, "Sm", j) * itensor::op(sites, "Id", j + 2) *
            itensor::op(sites, "Sp", j + 4);

      hh += -2 * J * itensor::op(sites, "Sp", j) * itensor::op(sites, "Sz", j + 2) *
            itensor::op(sites, "Sm", j + 4);

      hh += -2 * J * itensor::op(sites, "Sm", j) * itensor::op(sites, "Sz", j + 2) *
            itensor::op(sites, "Sp", j + 4);

      const double mu = j < dot ? hL * TL : hR * TR;

      hh += -mu * std::pow(-1, (j + 1) / 2) * itensor::op(sites, "Sz", j) *
            itensor::op(sites, "Id", j + 2) * itensor::op(sites, "Id", j + 4);
      if (j == end - 5) {
        hh += -mu * std::pow(-1, (j + 1 + 2) / 2) * itensor::op(sites, "Id", j) *
              itensor::op(sites, "Sz", j + 2) * itensor::op(sites, "Id", j + 4);
        hh += -mu * std::pow(-1, (j + 1 + 4) / 2) * itensor::op(sites, "Id", j) *
              itensor::op(sites, "Id", j + 2) * itensor::op(sites, "Sz", j + 4);
      }
      auto G = itensor::expHermitian(hh, -tau);
      gates.emplace_back(j, std::move(G));
    }
  }
  void TimeGates(const int begin, const int end, const std::complex<double> tau,
                 const itensor::SiteSet &sites, const ThreeSiteParam &param) {
    const int step = 6;
    const double J = param.value("J");
    for (int j = begin; j < end - 4; j += step) {
      std::cout << "j = (" << j << ", " << j + 2 << ", " << j + 4 << ")" << std::endl;
      // This part acts on physical sites.
      auto hh = J * itensor::op(sites, "Sp", j) * itensor::op(sites, "Id", j + 2) *
                itensor::op(sites, "Sm", j + 4);

      hh += J * itensor::op(sites, "Sm", j) * itensor::op(sites, "Id", j + 2) *
            itensor::op(sites, "Sp", j + 4);

      hh += -2 * J * itensor::op(sites, "Sp", j) * itensor::op(sites, "Sz", j + 2) *
            itensor::op(sites, "Sm", j + 4);

      hh += -2 * J * itensor::op(sites, "Sm", j) * itensor::op(sites, "Sz", j + 2) *
            itensor::op(sites, "Sp", j + 4);

      auto G = itensor::expHermitian(hh, -tau);
      gates.emplace_back(j, std::move(G));
    }
  }

  void EvolvePhysical(itensor::MPS &psi, const itensor::Args &args) {
    for (auto &gate : gates) {
      auto j = gate.i1;
      auto &G = gate.G;
      SwapNextSites(psi, j);     // swap j,j+1
      SwapNextSites(psi, j + 3); // Now physical sites are j+1,j+2,j+3
      psi.position(j + 2);
      auto WF = psi(j + 1) * psi(j + 2) * psi(j + 3);
      WF = G * WF;
      WF /= itensor::norm(WF);
      WF.noPrime();
      {
        auto [Uj1, Vj1] = itensor::factor(
            WF, {itensor::siteIndex(psi, j + 1), itensor::leftLinkIndex(psi, j + 1)}, args);
        auto indR = itensor::commonIndex(Uj1, Vj1);
        auto [Uj2, Vj2] = itensor::factor(Vj1, {itensor::siteIndex(psi, j + 2), indR}, args);
        psi.set(j + 1, Uj1);
        psi.set(j + 2, Uj2);
        psi.set(j + 3, Vj2);
        SwapNextSites(psi, j + 3);
        SwapNextSites(psi, j);
      }
    }
  }

  void Evolve(itensor::MPS &psi, const itensor::Args &args) { EvolvePhysical(psi, args); }
  void SwapNextSites(itensor::MPS &psi, const int j) {
    psi.position(j);
    auto WF = psi(j) * psi(j + 1);
    auto [U, V] = itensor::factor(
        WF, {itensor::siteIndex(psi, j + 1), itensor::leftLinkIndex(psi, j)}, {"Truncate=", false});
    psi.set(j, U);
    psi.set(j + 1, V);
    psi.position(j);
  }

private:
  std::vector<TGate> gates;
};

void DisconnectChains(itensor::MPS &psi, const int j) {
  psi.position(j);
  auto WF = psi(j) * psi(j + 1);
  auto [U, D, V] =
      itensor::svd(WF, {itensor::siteIndex(psi, j), itensor::leftLinkIndex(psi, j)}, {"MaxDim", 1});
  psi.set(j, U);
  psi.set(j + 1, D * V);
  psi.orthogonalize();
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
  itensor::MPS psi;
  std::cout << "Constructing Hamiltonian" << std::endl;
  ThreeSiteHamiltonian Ham(sites, param);
  const int dot = Ham.dot;
  auto H = itensor::toMPO(Ham.ampo);
  std::cout << "Finished constructing Hamiltonian" << std::endl;
  auto energy(0);
  psi = itensor::MPS(sites);
  std::cout << "N= " << N << std::endl;
  for (int i = 1; i < N; i += 2) { // N = 2L, where L is the original system size
    // Odd sites are physical sites and even sites are their ancillas.
    auto il = itensor::commonIndex(psi(i), psi(i + 1), "Link"); // Just to initialize the variable;
    auto ir = itensor::commonIndex(psi(i), psi(i + 1), "Link"); // Just to initialize the variable;
    auto s1 = sites(i);
    auto s2 = sites(i + 1);
    auto wf = itensor::ITensor(); // Just to initialize the variable;
    if (i == 1) {
      ir = itensor::commonIndex(psi(i + 1), psi(i + 2), "Link");
      wf = itensor::ITensor(s1, s2, ir);
      psi.ref(i) = itensor::ITensor(s1);
      psi.ref(i + 1) = itensor::ITensor(s2, ir);
      wf.set(s1(1), s2(2), ir(1), itensor::ISqrt2); // | +- > - | -+ > state
      wf.set(s1(2), s2(1), ir(1), -itensor::ISqrt2);
    } else if (i == N - 1) {
      il = itensor::commonIndex(psi(i - 1), psi(i), "Link");
      wf = itensor::ITensor(il, s1, s2);
      psi.ref(i) = itensor::ITensor(s1, il);
      psi.ref(i + 1) = itensor::ITensor(s2);
      wf.set(il(1), s1(1), s2(2), itensor::ISqrt2); // | +- > - | -+ > state
      wf.set(il(1), s1(2), s2(1), -itensor::ISqrt2);
    } else {
      il = itensor::commonIndex(psi(i - 1), psi(i), "Link");
      ir = itensor::commonIndex(psi(i + 1), psi(i + 2), "Link");
      wf = itensor::ITensor(il, s1, s2, ir);
      psi.ref(i) = itensor::ITensor(s1, il);
      psi.ref(i + 1) = itensor::ITensor(s2, ir);
      wf.set(il(1), s1(1), s2(2), ir(1), itensor::ISqrt2); // | +- > - | -+ > state
      wf.set(il(1), s1(2), s2(1), ir(1), -itensor::ISqrt2);
    }
    itensor::ITensor D;
    itensor::svd(wf, psi.ref(i), D, psi.ref(i + 1)); // put the prepared ancilla "wf" to the |PSI>
    psi.ref(i) *= D;
  }
  energy = itensor::inner(psi, H, psi); //<psi|H0|psi>
  std::cout << "2. Initial energy=" << energy << " .Norm = " << itensor::inner(psi, psi)
            << std::endl;

  double tau = param.value("tau");

  auto args = itensor::Args("Method=", "DensityMatrix", "Cutoff", param.value("trunc"), "MaxDim",
                            param.integer_value("max_bond"), "Normalize",
                            false); // arguments for time dynamics

  auto args0 = itensor::Args("Method=", "DensityMatrix", "Cutoff", param.value("trunc0"), "MaxDim",
                             param.integer_value("max_bond"), "Normalize",
                             false); // arguments for IMAGINARY time == initial temperature state

  // Output files for enabled observables.
  std::ofstream ent, spec, eprof, sz, sz_avrg, energy_beta, energy_prof, q1minus_prof,
      q2prof; // here I'm defining output streams == files
  std::ios_base::openmode mode;
  mode = std::ofstream::out; // Erase previous file (if present)

  double dt = param.value("Entropy");
  if (dt != 0) { // Entropy in the center of the chain
    ent.open("Entropy_center.dat", mode);
    ent.precision(15);
    ent << "#time \t Entropy(dot) \t BondDim(dot) \t MaxBondDim\n";
  }
  dt = param.value("SVD_spec");
  if (dt > 0) { // SVD Spectrum on central bond
    spec.open("SVD_spec.dat", mode);
    spec.precision(15);
    spec << "#Position=" << dot << "\t<SVD_spectrum>\t\ttime\n";
  }
  dt = param.value("Energy_beta");
  if (dt > 0) { // total Energy vs beta. plus std::max bond dim at the 3rd column
    energy_beta.open("Energy_beta.dat", mode);
    energy_beta.precision(15);
    energy_beta << "beta \t Energy \n";
  }
  dt = param.value("Eprof");
  if (dt > 0) { // Full entropy profile
    eprof.open("Entropy_profile.dat", mode);
    eprof.precision(15);
    eprof << "#Position=i-" << dot << std::setw(16) << "\t Entropy(i)" << std::setw(16)
          << "\t Entropy_sqrt \t Entropy_state1 \t time \t\t Bond.Dim(i)\n";
  }
  dt = param.value("Sz");
  if (dt > 0) { // Full magnetization profile
    sz.open("Sz_profile.dat", mode);
    sz.precision(15);
    sz << "#Position=i-" << "\t<Sz_i>\t" << "\t(-1)^i<Sz_i>\t" << "\t\ttime\n";

    sz_avrg.open("Sz_average_profile.dat", mode);
    sz_avrg.precision(15);
    sz_avrg << "#Position=i-" << "\t0.5<Sz_2+1i> + 0.5<Sz_2+1i>\t" << "\t\ttime\n";
  }
  dt = param.value("EnergyProf");
  if (dt > 0) { // Energy profile
    energy_prof.open("Energy_profile.dat", mode);
    energy_prof.precision(15);
    energy_prof << "#Position=i-" << "\t<Ham_i>\t" << dot << "\t\ttime(or beta)\n";

    // Q1minus profile is initialized with the energy profile.
    q1minus_prof.open("Q1minus_profile.dat", mode);
    q1minus_prof.precision(15);
    q1minus_prof << "#Position=i-" << "\t<Q1minus_i>\t" << dot << "\t\ttime(or beta)\n";
  }
  dt = param.value("Q2Prof");
  if (dt > 0) { // Full entropy profile
    q2prof.open("Q2_profile.dat", mode);
    q2prof.precision(15);
    q2prof << "#Position=i-" << dot << std::setw(16) << "\t Entropy(i)" << std::setw(16)
           << "\t Q2plus \t Q2minus \t time \t \n";
  }
  std::cout << "Trotter Gates for beta " << std::endl;
  param.set("begin", 1);
  param.set("end", N);
  const double dbeta = param.value("dbeta");
  TrotterExp expH_beta(sites, param, 0.5 * dbeta);

  std::cout << "Trotter Gates Half for beta " << std::endl;
  param.set("begin", dot + 1);
  TrotterExp expH_beta_half(sites, param, 0.5 * dbeta);

  std::cout << "Trotter Gates for tau" << std::endl;
  param.set("begin", 1);
  param.set("hL", 0);
  param.set("hR", 0);
  TrotterExp expH(sites, param, itensor::Cplx_i * 1.0 * tau);

  // Prepare the biased thermal state.
  const double TL = param.value("TL");
  const double TR = param.value("TR");
  const double beta_min = std::min(1. / TL, 1. / TR);
  const double beta_max = std::max(1. / TL, 1. / TR);
  const int beta_steps_min = beta_min / param.value("dbeta");
  const int beta_steps_max = beta_max / param.value("dbeta");
  const double n_steps = param.value("T") / param.value("tau");
  param.set("begin", 1);
  param.set("hL", 0);
  param.set("hR", 0);
  std::cout << "Finished preparing half-chain gates" << std::endl;

  const int time_total = beta_steps_max + n_steps;
  for (int n = 0; n <= time_total; ++n) {
    double time = (n - beta_steps_max) * tau; //+param.value("time_shift");

    if (n < beta_steps_max) {
      time = n * dbeta;
      std::cout << "Beta(1/T) step #" << n << "/" << time_total << "\t beta=" << time << std::endl;
    } else {
      std::cout << "Time step #" << n << "/" << time_total << "\ttime=" << time << std::endl;
    }
    std::cout.flush();
    std::vector<double> Myspec; // std::vector which will be the SVD spectrum

    if (param.value("Entropy") != 0) {
      double entr = Entropy(psi, dot, Myspec, 1);
      ent << time << "\t" << std::setw(16) << std::setfill('0') << entr << "\t" << BondDim(psi, dot)
          << "\t" << itensor::maxLinkDim(psi) << std::endl;

      if (param.value("SVD_spec") > 0) {
        if (n % int(param.value("SVD_spec") / tau) == 0) {
          spec << "\"t=" << time << "\"" << std::endl;
          int si = Myspec.size();
          for (int i = 0; i < si; i++) {
            spec << i + 1 << "\t" << Myspec[i] << "\t" << "\t" << time << std::endl;
          }
          if (n < time_total)
            spec << std::endl << std::endl;
        }
      }
    }

    if (param.value("Eprof") > 0)
      if (n % int(param.value("Eprof") / tau) == 0) {
        eprof << "\"t=" << time << "\"" << std::endl;
        for (int i = 1; i < N; i++) {
          double entr_std = Entropy(psi, i, Myspec, 1); // p log p
          eprof << i + 0.5 - dot << "\t" << std::setw(16) << std::setfill('0') << entr_std << "\t"
                << std::setw(4) << std::setfill('0') << BondDim(psi, i) << "\t" << time
                << std::endl;
        }
        if (n < time_total)
          eprof << "\n\n"; // I need this part to separate time steps in gnuplot
      }

    if (param.value("Energy_beta") > 0) {
      double en = 0;
      double counter = 0;
      int shift = (N / 4) % 2 == 0 ? N / 4 : N / 4 + 1;
      for (int i = dot + 1 - shift; i < dot + 1 + shift; i += 2) {
        en += Energy(psi, sites, i);
        ++counter;
      }
      en /= (counter);

      energy_beta << time << "\t"
                  << std::real(itensor::innerC(psi, H, psi)) / std::real(itensor::innerC(psi, psi))
                  << "\t" << en << "\t" << itensor::maxLinkDim(psi) << std::endl;
    }
    if (param.value("EnergyProf") > 0 && beta_steps_max <= n) {
      if (n % int(param.value("EnergyProf") / tau) == 0) {
        energy_prof << "\"t=" << time << "\"" << std::endl;
        q1minus_prof << "\"t=" << time << "\"" << std::endl;
        for (int i = 1; i <= N - 5; i += 2) {
          const std::complex<double> q1 = Q1(psi, sites, i);
          const double en = std::real(q1);
          energy_prof << i / 2 - dot / 2 + 1 << "\t" << en << "\t" << time << std::endl;
          const double q1minus = std::imag(q1);
          q1minus_prof << i / 2 - dot / 2 + 1 << "\t" << q1minus << "\t" << time << std::endl;
        }
        energy_prof
            << "\n\n"; // I need this part to separate time steps in *.dat files (for gnuplot)
        q1minus_prof
            << "\n\n"; // I need this part to separate time steps in *.dat files (for gnuplot)
      }
    }

    if (param.value("Q2Prof") > 0 && beta_steps_max <= n) {
      if (n % int(param.value("EnergyProf") / tau) == 0) {
        q2prof << "\"t=" << time << "\"" << std::endl;
        for (int i = 1; i <= N - 9; i += 2) {
          const std::complex<double> q2 = Q2(psi, sites, i);
          const double q2plus = std::real(q2);
          const double q2minus = std::imag(q2);
          q2prof << i / 2 - dot / 2 + 1 << "\t" << q2plus << "\t" << q2minus << "\t" << time
                 << std::endl;
        }
        q2prof << "\n\n"; // I need this part to separate time steps in *.dat files (for gnuplot)
      }
    }

    if (param.value("Sz") > 0 && beta_steps_max <= n) {
      if (n % int(param.value("Sz") / tau) == 0) {
        sz << "\"t=" << time << "\"" << std::endl;
        double sz_tot = 0, sz_left = 0, sz_right = 0, sz_dot = 0;
        double sz_odd = 0;
        for (int i = 1; i <= N; i += 2) {
          double s = Sz(psi, sites, i);
          sz_tot += s;
          if (i < dot)
            sz_left += s;
          if (i > dot)
            sz_right += s;
          if (i == dot)
            sz_dot += s;
          sz << i / 2 - dot / 2 + 1 << "\t" << s << "\t" << std::pow(-1, (i + 1) / 2) * s << "\t"
             << time << std::endl;

          if (((i + 1) / 2) % 2 == 1) { // odd site
            sz_odd = s;
          } else {
            sz_avrg << i / 2 - dot / 2 + 1 << "\t" << 0.5 * (s + sz_odd) << "\t" << time
                    << std::endl;
          }
        }

        { // I need this part to separate time steps in *.dat files (for gnuplot)
          sz << "\n\n";
          sz_avrg << "\n\n";
        }
        std::cout << "\n<Sz_left>=" << sz_left << "\t" << "<Sz_right>=" << sz_right << "\t"
                  << "<Sz_DOT>=" << sz_dot << "\t" << "<Sz_tot>=" << sz_tot << std::endl;
      }
    }

    if (n < beta_steps_min) {
      std::cout << "Temperature evolution with H" << std::endl;
      expH_beta.EvolvePhysical(psi, args0);
      psi.orthogonalize(args);
      std::cout << "dot = " << dot + 1 << std::endl;

      psi.normalize();
    } else if (n < beta_steps_max) {
      std::cout << "Temperature evolution with half-chain H" << std::endl;
      std::cout << "n/beta_steps_max = " << n << "/" << beta_steps_max << std::endl;
      expH_beta_half.EvolvePhysical(psi, args0);
      psi.orthogonalize(args);
      psi.normalize();
    } else {
      std::cout << "Time evolution" << std::endl;
      expH.Evolve(psi, args);
      psi.orthogonalize(args);
    }
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
