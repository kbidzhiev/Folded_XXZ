#pragma once

#include <iosfwd>
#include <map>
#include <string>

class Parameters {
public:
  double value(const std::string &name) const;
  int integer_value(const std::string &name) const;
  void set(const std::string &name, double value);
  bool contains(const std::string &name) const;
  void print(std::ostream &os) const;
  void parse_arguments(int argc, char *argv[]);

protected:
  void define(const std::string &name, double value);

private:
  std::map<std::string, double> values_;
};

class ThreeSiteParam : public Parameters {
public:
  ThreeSiteParam();
};
