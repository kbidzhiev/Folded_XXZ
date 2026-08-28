#include <folded_xxz/parameters.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {

double parse_double(char *text) {
  char *end = nullptr;
  const double value = std::strtod(text, &end);
  if (end == text || *end != '\0') {
    std::cerr << "ERROR: " << text << " is not a valid double." << std::endl;
    std::exit(1);
  }
  return value;
}

} // namespace

double Parameters::value(const std::string &name) const {
  const auto it = values_.find(name);
  if (it == values_.end())
    throw std::out_of_range("unknown parameter: " + name);
  return it->second;
}

int Parameters::integer_value(const std::string &name) const {
  const double number = value(name);
  const double rounded = std::round(number);
  if (std::abs(rounded - number) > 1e-6)
    throw std::invalid_argument("parameter " + name + " must be an integer");
  if (rounded < std::numeric_limits<int>::min() || rounded > std::numeric_limits<int>::max())
    throw std::out_of_range("parameter " + name + " is outside the int range");
  return static_cast<int>(rounded);
}

void Parameters::set(const std::string &name, double value) {
  if (!contains(name))
    throw std::out_of_range("unknown parameter: " + name);
  values_.at(name) = value;
}

bool Parameters::contains(const std::string &name) const {
  return values_.find(name) != values_.end();
}

void Parameters::print(std::ostream &os) const {
  for (const auto &[name, value] : values_)
    os << name << "=" << value << std::endl;
}

void Parameters::parse_arguments(int argc, char *argv[]) {
  for (int index = 1; index < argc; ++index) {
    const std::string name(argv[index]);
    if (!contains(name)) {
      std::cerr << "Syntax error: " << name << std::endl;
      print(std::cerr);
      std::exit(1);
    }
    if (++index == argc) {
      std::cerr << "Error: missing value after " << name << std::endl;
      std::exit(1);
    }
    set(name, parse_double(argv[index]));
  }
}

void Parameters::define(const std::string &name, double value) { values_.emplace(name, value); }

ThreeSiteParam::ThreeSiteParam() {
  define("N", 10);
  define("J", 1.0);
  define("tau", 0.01);
  define("dbeta", 0.01);
  define("T", 0);
  define("TL", 100);
  define("TR", 5);
  define("hL", 0);
  define("hR", 0);
  define("Entropy", 0);
  define("Eprof", 0);
  define("EnergyProf", 0);
  define("Q2Prof", 0);
  define("Sz", 0);
  define("SVD_spec", 0);
  define("max_bond", 4000);
  define("trunc", 1e-10);
  define("trunc0", 1e-10);
  define("TrotterOrder", 2);
  define("Energy_beta", 1);
}
