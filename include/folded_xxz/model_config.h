#pragma once

struct ModelConfig {
  int site_count;
  double coupling;
  double left_thermal_field;
  double right_thermal_field;
  int trotter_order;
};
