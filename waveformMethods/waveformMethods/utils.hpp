#ifndef WAVEFOMR_METHODS_UTILS_H
#define WAVEFOMR_METHODS_UTILS_H

#include "waveformMethods/waveformMethods.hpp"

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace waveform_methods::Utils {

  double LinearInterpolationX(
    const double &x1, const double &y1,
    const double &x2, const double &y2,
    const double &y);

  template<class dtype>
  dtype MedienFilter(
    const std::vector<dtype> &data,
    const int &window_size);

}

#endif
